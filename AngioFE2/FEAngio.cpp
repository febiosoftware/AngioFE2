///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FECore/FEDataLoadCurve.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMix/FESolute.h"
#include "FEBioMix/FEMultiphasic.h"
#include "FEBioLib/FEBioModel.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FECore/FESolidDomain.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "Elem.h"
#include "angio3d.h"
#include <ctime>
#include <future>
#include <algorithm>
#include <cfloat> 
#include <FEBioMix/FETriphasic.h>


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat);

bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID);

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : ztopi(std::uniform_real_distribution<double>(0, pi)), 
	zto2pi(std::uniform_real_distribution<double>(0, 2 * pi)), n1to1(std::uniform_real_distribution<double>(-1, 1))
{
	// Body force counter
	total_bdyf = 0;
	
	FE_state = 0;

	// initialize time stepping parameters
	m_time.dt = 1.0;
	// Input random seed number
	m_irseed = 0;

	m_fem = dynamic_cast<FEBioModel *>(&fem);
	assert(m_fem);
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
}

//-----------------------------------------------------------------------------
FEBioModel* FEAngio::GetFEModel() const
{
	return m_fem;
}

FEMesh * FEAngio::GetMesh() const
{
	return &m_fem->GetMesh();
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	//create any classes which have nontrivial destructors 
	//currently the destructors are not called for classes created by FEBio this allows destructors to be called
	fileout = new Fileout();

	// Init all the FE stuff
	//must be done first initializes material
	if (InitFEM() == false) return false;

	// Seed the random number generator based on the sum of the seeds of the materials ie set the seed only once in any material
	unsigned int posseed = static_cast<unsigned int>(m_fem->GetGlobalConstant("seed"));
	if (!posseed)
	{
		posseed = static_cast<unsigned int>(time(0));
	}
	m_irseed = posseed;
	rengine.seed(m_irseed);
	srand(m_irseed);

	ForEachNode([this](FENode & node)
	{
		//nodes may be accessed multiple times by the current implementation
		FEAngioNodeData nd;
		m_fe_node_data[node.GetID()] = nd;
	});
	ForEachElement([this](FEElement & e, FEDomain & d)
	               {
		               if (m_fe_element_data.count(e.GetID())) assert(false);
		               FEAngioElementData ed;
		               m_fe_element_data[e.GetID()] = ed;
	               });


	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->FinalizeInit();
	}
	std::vector<std::future<void>> surface_futures;
	//setup the exterior_surface
	SetupSurface(surface_futures);

	// assign ECM densities to grid nodes
	//TODO: degree of anisotropy values will only get initialized if density is set to 0
	if (InitECMDensity() == false) return false;

	CallInFutures(surface_futures);
	// assign concentration values to grid nodes
	//if (InitSoluteConcentration() == false) return false;

	// NOTE: must be done after InitECMDensity() and InitCollagenFibers().
	bool rv = true;
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		rv &= m_pmat[i]->InitCulture();
	}
	if (!rv)
		return false;

	FinalizeFEM();
	


	// start timer
	time(&m_start);

	return true;
}

void FEAngio::SetupSurface(std::vector<future<void>> & futures)
{
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		futures.push_back(std::async(std::launch::async,&FEAngioMaterialBase::SetupSurface, m_pmat[i]));
	}
}
void FEAngio::CallInFutures(std::vector<future<void>> & futures)
{
	for(int i =0; i < futures.size();i++)
	{
		futures[i].get();
	}
}

//-----------------------------------------------------------------------------
double FEAngio::RunTime() const
{
	time_t stop;
	time(&stop);
	return static_cast<double>(difftime(stop, m_start));
}

//-----------------------------------------------------------------------------
// Initialize FE model.
bool FEAngio::InitFEM()
{
	for (int i = 0; i < m_fem->Materials(); i++)
	{
		FEAngioMaterial * cmat = nullptr;
		FEMaterial * mat = m_fem->GetMaterial(i);
		int id = -1;
		cmat = GetAngioComponent(mat);
		

		if (cmat)
		{
			id = mat->GetID();
			assert(id != -1);
			m_pmat.emplace_back(cmat);
			m_pmat_ids.emplace_back(id);
			//TODO: check that material parameters are set here
			//cmat->ApplySym();
			cmat->SetFEAngio(this);
			cmat->UpdateSproutStressScaling();
		}
	}
	//assert(m_pmat.size());

	felog.printf("%d Angio materials found. Stress approach will be used.", m_pmat.size());

	// register the callback
	m_fem->AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED | CB_STEP_ACTIVE, this);

	// Do the model initialization
	if (m_fem->Init() == false) return false;

	return true;
}
void FEAngio::FinalizeFEM()
{
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		//TODO: remove this constant or make it a user parameter
		m_pmat[i]->CreateSprouts(0.5,m_pmat[i]->GetMatrixMaterial()->GetElasticMaterial());
		m_pmat[i]->AdjustMeshStiffness(m_pmat[i]->GetMaterial());
		m_pmat[i]->UpdateFiberManager();
		m_pmat[i]->InitializeFibers();
		m_pmat[i]->UpdateFiberManager();
	}

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::LOG_FILE);


	// --- Output initial state of model ---
	if (!m_fem->GetGlobalConstant("no_io"))
	{
		// Output initial microvessel state
		fileout->save_vessel_state(*this);

		// save active tips
		fileout->save_active_tips(*this);
	}
}

//-----------------------------------------------------------------------------
// Initialize the nodal ECM values
bool FEAngio::InitECMDensity()
{
	ForEachNodePar([&](FENode & node)
	{
		m_fe_node_data[node.GetID()].m_collfib = vec3d(0, 0, 0);
		m_fe_node_data[node.GetID()].m_ecm_den = 0.0;
	});
	bool rv = true;
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		rv &= m_pmat[i]->InitECMDensity(this);
	}
	if (!rv)
	{
		return false;
	}

	// normalize fiber vector and average ecm density
	ForEachNodePar([&](FENode & node)
	{
		//nneds to be run only once per node
		if (m_fe_node_data[node.GetID()].m_ntag)
		{
			m_fe_node_data[node.GetID()].m_ecm_den = m_fe_node_data[node.GetID()].m_ecm_den0;
			m_fe_node_data[node.GetID()].m_ntag = 0;
		}
	});
	return true;
}

double FEAngio::GetDoubleFromDataStore(int record, int elem_id, int item)
{
	DataStore & ds = m_fem->GetDataStore();
	return ds.GetDataRecord(record)->Evaluate(elem_id, item);
}

//-----------------------------------------------------------------------------
// update the extracellular matrix nodal values
void FEAngio::UpdateECM()
{
	//straight translation of code from the grid probably should be optimized later

	// reset nodal data
	FEMesh & mesh = m_fem->GetMesh();

	ForEachNodePar([&](FENode & node)
	{
		m_fe_node_data[node.GetID()].m_collfib = vec3d(0, 0, 0);
		m_fe_node_data[node.GetID()].m_ecm_den = 0.0;
		m_fe_node_data[node.GetID()].m_ntag = 0;
		//REFACTOR: why reset not just overwrite
	});

	//this portion will be harder
	// For each element within the grid...
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->UpdateECM();
	}

	// normalize fiber vector and average ecm density
	ForEachNode([this](FENode & node)
	{
		//nneds to be run only once per node
		if (m_fe_node_data[node.GetID()].m_ntag)
		{
			m_fe_node_data[node.GetID()].m_ecm_den /= (double)m_fe_node_data[node.GetID()].m_ntag;
			m_fe_node_data[node.GetID()].m_collfib.unit();
			m_fe_node_data[node.GetID()].m_ntag = 0;
			//maybe worry about density creep
			//m_fe_node_data[node.GetID_ang()].m_ecm_den0 = m_fe_node_data[node.GetID_ang()].m_ecm_den;
		}
	});
}

int FEAngio::FindGrowTimes(std::vector<std::pair<double, double>> & time_pairs, int start_index)
{
	double res = m_fem->GetGlobalConstant("angio_time_step_curve");
	//0 will be returned on failure
	int curve_index = static_cast<int>(res);
	if (curve_index)
	{
		FEDataLoadCurve * fecurve = dynamic_cast<FEDataLoadCurve*>(m_fem->GetLoadCurve(curve_index -1));
		SimulationTime st = CurrentSimTime();
		assert(fecurve);
		int index = start_index;
		double start_time = st.t - st.dt;
		for (int i = start_index; i < fecurve->Points(); i++)
		{
			LOADPOINT lp = fecurve->LoadPoint(i);
			if (lp.time <= st.t)
			{
				double dt = lp.time - start_time;
				assert(dt > 0.0);
				time_pairs.emplace_back(start_time, dt);
				start_time = lp.time;
			}
			else
			{
				return i;
			}
		}

		return index;
	}
	else
	{
		SimulationTime st = CurrentSimTime();
		time_pairs.emplace_back(st.t-st.dt, st.dt);
	}
	return 0;
}


//TODO: consider making the ForEach* const as long as this works on both compilers
void FEAngio::ForEachNode(std::function<void(FENode &)> f, std::vector<int> & matls)
{
	//TODO: the last element to access a node wins on overwting the data associated with that node
	//this behavior matches the previous behavior of the plugin but probably should be fixed sometime
	FEMesh & mesh = m_fem->GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	assert(matls.size() == 0 || dl.size());
	for (size_t i = 0; i < dl.size(); i++)
	{
		FEDomain & d = mesh.Domain(dl[i]);
		for (int j = 0; j < d.Elements(); j++)
		{
			FEElement & e = d.ElementRef(j);
//#pragma omp parallel for schedule(dynamic)
			for (int k = 0; k < e.Nodes(); k++)
			{
				f(mesh.Node(e.m_node[k]));//this iterates over the local nodes
			}
		}
	}
}

void FEAngio::ForEachNodePar(std::function<void(FENode &)> f, std::vector<int> & matls)
{
	//TODO: the last element to access a node wins on overwting the data associated with that node
	//this behavior matches the previous behavior of the plugin but probably should be fixed sometime
	FEMesh & mesh = m_fem->GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	assert((matls.size() ==0 )|| dl.size());
	for (size_t i = 0; i < dl.size(); i++)
	{
		FEDomain & d = mesh.Domain(dl[i]);
		for (int j = 0; j < d.Elements(); j++)
		{
			FEElement & e = d.ElementRef(j);
#pragma omp parallel for schedule(dynamic)
			for (int k = 0; k < e.Nodes(); k++)
			{
				f(mesh.Node(e.m_node[k]));//this iterates over the local nodes
			}
		}
	}
}
void FEAngio::ForEachNodePar(std::function<void(FENode &)> f)
{
	ForEachNodePar(f, m_pmat_ids);
}
void FEAngio::ForEachNode(std::function<void(FENode &)> f)
{
	ForEachNode(f, m_pmat_ids);
}
void FEAngio::ForEachElement(std::function<void(FESolidElement&, FESolidDomain&)> f, std::vector<int> & matls)
{
	FEMesh & mesh = m_fem->GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	assert(matls.size() == 0 || dl.size());
	for (size_t i = 0; i < dl.size(); i++)
	{
		FESolidDomain & d = reinterpret_cast<FESolidDomain&>(mesh.Domain(dl[i]));
		for (int j = 0; j < d.Elements(); j++)
		{
			FESolidElement & e = reinterpret_cast<FESolidElement&>(d.ElementRef(j));
			f(e, d);
		}
	}
}
void FEAngio::ForEachElementPar(std::function<void(FESolidElement&, FESolidDomain&)> f, std::vector<int> & matls)
{
	FEMesh & mesh = m_fem->GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	assert(matls.size() ==0 || dl.size());
	for (size_t i = 0; i < dl.size(); i++)
	{
		FESolidDomain & d = reinterpret_cast<FESolidDomain&>(mesh.Domain(dl[i]));
#pragma omp parallel for schedule(dynamic, 24)
		for (int j = 0; j < d.Elements(); j++)
		{
			FESolidElement & e = reinterpret_cast<FESolidElement&>(d.ElementRef(j));
			f(e, d);
		}
	}
}
void FEAngio::ForEachElementPar(std::function<void(FESolidElement&, FESolidDomain&)> f)
{
	ForEachElementPar(f, m_pmat_ids);
}
void FEAngio::ForEachElement(std::function<void(FESolidElement&, FESolidDomain&)> f)
{
	ForEachElement(f, m_pmat_ids);
}
void FEAngio::ForEachDomain(std::function<void(FESolidDomain&)> f)
{
	ForEachDomain(f, m_pmat_ids);
}
void FEAngio::ForEachDomain(std::function<void(FESolidDomain&)> f, std::vector<int> & matls)
{
	FEMesh & mesh = m_fem->GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	for (size_t i = 0; i < dl.size(); i++)
	{
		FESolidDomain & d = reinterpret_cast<FESolidDomain&>(mesh.Domain(dl[i]));
		f(d);
	}
}

mat3d FEAngio::unifromRandomRotationMatrix()
{
	//collagen fibers are right handed so the following transformation is legal
	//the following will only produce right handed bases for the collagen fibers which is molecularly accurate
	double alpha = zto2pi(rengine);
	double beta = zto2pi(rengine);
	double gamma = zto2pi(rengine);


	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, - s_alpha*s_gamma+c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
		);
	//rv = rv.transinv();
#ifndef NDEBUG
	//do some testing that bases run tthrough this are still orthogonal
	vec3d x(1, 0, 0);
	vec3d y(0, 1, 0);
	vec3d z(0, 0, 1);
	vec3d xt = rv * x;
	vec3d yt = rv * y;
	vec3d zt = rv * z;
	double tol = 0.01;
	assert(xt * yt < tol && xt * yt > -tol);
	assert(xt * zt < tol && xt * zt > -tol);
	assert(zt * yt < tol && zt * yt > -tol);

	double noq = xt.norm();
	double nom1 = yt.norm();
	double nom2 = zt.norm();
	tol = 0.01;
	assert(noq < (1 + tol) && noq >(1 - tol));
	assert(nom1 < (1 + tol) && nom1 >(1 - tol));
	assert(nom2 < (1 + tol) && nom2 >(1 - tol));
#endif
	return rv;
}
mat3d FEAngio::rotationMatrix(double alpha, double beta, double gamma)
{
	double c_alpha = cos(alpha);
	double c_beta = cos(beta);
	double c_gamma = cos(gamma);
	double s_alpha = sin(alpha);
	double s_beta = sin(beta);
	double s_gamma = sin(gamma);
	//see: https://en.wikipedia.org/wiki/Change_of_basis Three dimensions section
	mat3d rv(c_alpha*c_gamma - s_alpha*c_beta*s_gamma, -c_alpha*s_gamma - s_alpha*c_beta*c_gamma, s_beta*s_alpha,
		s_alpha*c_gamma + c_alpha*c_beta*s_gamma, -s_alpha*s_gamma + c_alpha*c_beta*c_gamma, -s_beta*c_alpha,
		s_beta*s_gamma, s_beta*c_gamma, c_beta
	);
	return rv;
}

vec3d FEAngio::uniformRandomDirection()
{
	//to revert this set this to return vrand
	double theta = zto2pi(rengine);
	double phi = zto2pi(rengine);
	double sintheta = sin(theta);
	vec3d dir(sintheta * cos(phi), sintheta * sin(phi), cos(theta));
	return dir;
}
vec3d FEAngio::uniformInUnitCube()
{
	//force the order of initialization should allow results on different platforms to converge 
	double x = n1to1(rengine);
	double y = n1to1(rengine);
	double z = n1to1(rengine);
	return vec3d(x,y,z);
}

vec3d FEAngio::gradient(FESolidElement * se, std::vector<double> & fn, vec3d pt, int size, int offset)
{
	assert(se);
	FESolidDomain* domain = dynamic_cast<FESolidDomain*>(se->GetDomain());
	double Ji[3][3];
	domain->invjact(*se, Ji, pt.x, pt.y, pt.z);
	double Gr[FEElement::MAX_NODES], Gs[FEElement::MAX_NODES], Gt[FEElement::MAX_NODES];
	se->shape_deriv(Gr, Gs, Gt, pt.x, pt.y, pt.z);

	vec3d gradf;
	size_t N = se->Nodes();
	assert(N == size*fn.size());
	for (size_t i = 0; i<N; ++i)
	{
		double Gx, Gy, Gz;
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0] * Gr[i] + Ji[1][0] * Gs[i] + Ji[2][0] * Gt[i];
		Gy = Ji[0][1] * Gr[i] + Ji[1][1] * Gs[i] + Ji[2][1] * Gt[i];
		Gz = Ji[0][2] * Gr[i] + Ji[1][2] * Gs[i] + Ji[2][2] * Gt[i];

		// calculate pressure gradient
		gradf.x += Gx*fn[i*size + offset];
		gradf.y += Gy*fn[i*size + offset];
		gradf.z += Gz*fn[i*size + offset];
	}

	return gradf;
}
//-----------------------------------------------------------------------------
static void solve_3x3(double A[3][3], double b[3], double x[3])
{
	double D = A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[2][1] * A[0][2] \
		- A[1][1] * A[2][0] * A[0][2] - A[2][2] * A[1][0] * A[0][1] - A[0][0] * A[2][1] * A[1][2];

	assert(D != 0);

	double Ai[3][3];
	Ai[0][0] = A[1][1] * A[2][2] - A[2][1] * A[1][2];
	Ai[0][1] = A[2][1] * A[0][2] - A[0][1] * A[2][2];
	Ai[0][2] = A[0][1] * A[1][2] - A[1][1] * A[0][2];

	Ai[1][0] = A[2][0] * A[1][2] - A[1][0] * A[2][2];
	Ai[1][1] = A[0][0] * A[2][2] - A[2][0] * A[0][2];
	Ai[1][2] = A[1][0] * A[0][2] - A[0][0] * A[1][2];

	Ai[2][0] = A[1][0] * A[2][1] - A[2][0] * A[1][1];
	Ai[2][1] = A[2][0] * A[0][1] - A[0][0] * A[2][1];
	Ai[2][2] = A[0][0] * A[1][1] - A[0][1] * A[1][0];

	x[0] = (Ai[0][0] * b[0] + Ai[0][1] * b[1] + Ai[0][2] * b[2]) / D;
	x[1] = (Ai[1][0] * b[0] + Ai[1][1] * b[1] + Ai[1][2] * b[2]) / D;
	x[2] = (Ai[2][0] * b[0] + Ai[2][1] * b[1] + Ai[2][2] * b[2]) / D;


#ifdef _DEBUG
	double r[3];
	r[0] = b[0] - (A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2]);
	r[1] = b[1] - (A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2]);
	r[2] = b[2] - (A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2]);

	double nr = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
#endif
}

bool FEAngio::IsInsideHex8(FESolidElement * se, vec3d y, FEMesh * mesh, double r[3])
{
	assert(se->Type() == FE_HEX8G8);
	vec3d x[FEElement::MAX_NODES];
	int j;
	// get the element nodal coordinates
	int neln = se->Nodes();
	for (j = 0; j<neln; ++j) x[j] = mesh->Node(se->m_node[j]).m_rt;

	// first, as a quick check, we see if y lies in the bounding box defined by x
	FEBoundingBox box(x[0]);
	for (j = 1; j<neln; ++j) box.add(x[j]);

	if (box.IsInside(y))
	{
		// If the point y lies inside the box, we apply a Newton method to find
		// the isoparametric coordinates r
		r[0] = r[1] = r[2] = 0;
		const double tol = 1e-5;
		double dr[3], norm;
		double H[8], G[8][3];
		do
		{
			H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
			H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
			H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
			H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
			H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
			H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
			H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
			H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

			G[0][0] = -0.125*(1 - r[1])*(1 - r[2]); G[0][1] = -0.125*(1 - r[0])*(1 - r[2]); G[0][2] = -0.125*(1 - r[0])*(1 - r[1]);
			G[1][0] = 0.125*(1 - r[1])*(1 - r[2]); G[1][1] = -0.125*(1 + r[0])*(1 - r[2]); G[1][2] = -0.125*(1 + r[0])*(1 - r[1]);
			G[2][0] = 0.125*(1 + r[1])*(1 - r[2]); G[2][1] = 0.125*(1 + r[0])*(1 - r[2]); G[2][2] = -0.125*(1 + r[0])*(1 + r[1]);
			G[3][0] = -0.125*(1 + r[1])*(1 - r[2]); G[3][1] = 0.125*(1 - r[0])*(1 - r[2]); G[3][2] = -0.125*(1 - r[0])*(1 + r[1]);
			G[4][0] = -0.125*(1 - r[1])*(1 + r[2]); G[4][1] = -0.125*(1 - r[0])*(1 + r[2]); G[4][2] = 0.125*(1 - r[0])*(1 - r[1]);
			G[5][0] = 0.125*(1 - r[1])*(1 + r[2]); G[5][1] = -0.125*(1 + r[0])*(1 + r[2]); G[5][2] = 0.125*(1 + r[0])*(1 - r[1]);
			G[6][0] = 0.125*(1 + r[1])*(1 + r[2]); G[6][1] = 0.125*(1 + r[0])*(1 + r[2]); G[6][2] = 0.125*(1 + r[0])*(1 + r[1]);
			G[7][0] = -0.125*(1 + r[1])*(1 + r[2]); G[7][1] = 0.125*(1 - r[0])*(1 + r[2]); G[7][2] = 0.125*(1 - r[0])*(1 + r[1]);

			double R[3] = { 0 }, A[3][3] = { 0 };
			for (j = 0; j<8; ++j)
			{
				R[0] += x[j].x*H[j];
				R[1] += x[j].y*H[j];
				R[2] += x[j].z*H[j];

				A[0][0] -= x[j].x*G[j][0]; A[0][1] -= x[j].x*G[j][1]; A[0][2] -= x[j].x*G[j][2];
				A[1][0] -= x[j].y*G[j][0]; A[1][1] -= x[j].y*G[j][1]; A[1][2] -= x[j].y*G[j][2];
				A[2][0] -= x[j].z*G[j][0]; A[2][1] -= x[j].z*G[j][1]; A[2][2] -= x[j].z*G[j][2];
			}
			R[0] = y.x - R[0];
			R[1] = y.y - R[1];
			R[2] = y.z - R[2];

			solve_3x3(A, R, dr);
			r[0] -= dr[0];
			r[1] -= dr[1];
			r[2] -= dr[2];

			norm = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
		} while (norm > tol);

		// see if the point r lies inside the element
		const double eps = 1.0001;
		if ((r[0] >= -eps) && (r[0] <= eps) &&
			(r[1] >= -eps) && (r[1] <= eps) &&
			(r[2] >= -eps) && (r[2] <= eps)) return true;
	}
	return false;
}

vec3d FEAngio::LocalToGlobal(FESolidElement * se, vec3d & rst) const
{
	vec3d cur = vec3d(0, 0, 0);
	FEMesh * mesh = GetMesh();
	double sf[FEElement::MAX_NODES];
	se->shape_fnc(sf,rst.x, rst.y, rst.z);
	for (int k = 0; k < se->Nodes(); k++)
	{
		cur += mesh->Node(se->m_node[k]).m_rt * sf[k];
	}
	return cur;
}

vec3d FEAngio::FindRST(const vec3d & v ,vec2d rs, FESolidElement * elem) const
{
	//try brute force it may yeild interesting results: ie the distance between the actual point on surface and the gridpoint
	std::vector<vec3d> positions_global;
	std::vector<vec3d> positions_local;
	vec3d cur;
	//r
	vec3d cpos = vec3d(1, rs[0], rs[1]);
	cur = LocalToGlobal(elem,cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(1, rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//s
	cpos = vec3d(rs[0], 1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], 1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[0], -1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], -1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//t
	cpos = vec3d(rs[0], rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[0], rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//do all of the negations
	//first
	//r
	cpos = vec3d(1, -rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(1, -rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, -rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, -rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//s
	cpos = vec3d(-rs[0], 1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], 1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[0], -1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], -1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//t
	cpos = vec3d(-rs[0], rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[0], rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//second
	//r
	cpos = vec3d(1, rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(1, rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//s
	cpos = vec3d(rs[0], 1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], 1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[0], -1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], -1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//t
	cpos = vec3d(rs[0], -rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], -rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[0], -rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(rs[1], -rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);


	//both
	//r
	cpos = vec3d(1, -rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(1, -rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, -rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-1, -rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//s
	cpos = vec3d(-rs[0], 1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], 1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[0], -1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], -1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	//t
	cpos = vec3d(-rs[0], -rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], -rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[0], -rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);

	cpos = vec3d(-rs[1], -rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.emplace_back(cur);
	positions_local.emplace_back(cpos);


	double best_dist = DBL_MAX;
	size_t best_index = -1;
	for (size_t i = 0; i < positions_global.size(); i++)
	{
		double cdist = (v - positions_global[i]).norm();
		if (cdist < best_dist)
		{
			best_index = i;
			best_dist = cdist;
		}
	}
	if (best_dist > 0.25)
	{
		//may need negations if the method does not appear to work
		printf("large error here\n");
		assert(false);
	}
	return positions_local[best_index];
}

GridPoint FEAngio::FindGridPoint(FESolidDomain * domain, int nelem, vec3d& q) const
{
	assert(domain != nullptr && nelem >= 0);
	GridPoint pt;
	pt.q = q;
	FESolidElement * se;
	//TODO: refactor if problems with multiple domains
	if (se = &domain->Element(nelem))
	{
		pt.ndomain = domain;//set the domain of the gridpoint
		pt.elemindex = nelem;
		pt.nelem = se->GetID();
	}
	else
	{
		assert(false);
	}
	pt.r = Position(pt);
	return pt;
}

vec3d FEAngio::Position(const GridPoint& pt) const
{
	//Point has already been positioned
	FEMesh & mesh = m_fem->GetMesh();
	vec3d r(0, 0, 0);
	FESolidDomain * d = pt.ndomain;
	FESolidElement * se;
	if (se = &d->Element(pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, pt.q.x, pt.q.y, pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			r += mesh.Node(se->m_node[j]).m_rt* arr[j];
		}
	}
	return r;
}

vec3d FEAngio::ReferenceCoords(const GridPoint& pt) const
{
	vec3d r(0, 0, 0);
	FEMesh & mesh = m_fem->GetMesh();
	//Point has already been positioned
	FESolidElement * se;
	if (se = &pt.ndomain->Element(pt.elemindex))
	{
		double arr[FEElement::MAX_NODES];
		se->shape_fnc(arr, pt.q.x, pt.q.y, pt.q.z);
		for (int j = 0; j < se->Nodes(); j++)
		{
			r += mesh.Node(se->m_node[j]).m_r0* arr[j];
		}
	}
	return r;
}
FEAngioMaterial * FEAngio::GetAngioComponent(FEMaterial * mat)
{
	auto angm = dynamic_cast<FEAngioMaterial*>(mat);
	if (!angm)
	{
		FEMultiphasic * mmat = dynamic_cast<FEMultiphasic*>(mat);
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	if (!angm)
	{
		FETriphasic * mmat = dynamic_cast<FETriphasic*>(mat);
		if(mmat)
			angm = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
	}
	/* //Biphasic Material doesn't have a get solid function 
	 * //this fucntion should be moved to an interface in febio
	if (!angm)
	{
		FEBiphasic * mmat = dynamic_cast<FEBiphasic*>(mat);
		angm = dynamic_cast<FEAngioMaterialBase*>(mmat->GetSolid());
	}
	*/
	return angm;
}

//TODO: consider having multiple projection methods
std::vector<double> FEAngio::createVectorOfMaterialParameters(FEElement * elem,
	double FEAngioNodeData::*materialparam)
{
	FEMesh * mesh = GetMesh();
	std::vector<double> gx(elem->m_node.size());
	for (size_t i = 0; i < elem->m_node.size(); i++)
	{
		gx[i] = this->m_fe_node_data[mesh->Node(elem->m_node[i]).GetID()].*materialparam;
	}
	return gx;
}
double FEAngio::genericProjectToPoint(FESolidElement * elem,
	double FEAngioNodeData::*materialparam,const vec3d & pos)
{
	assert(elem);
	std::vector<double> gx = createVectorOfMaterialParameters( elem, materialparam);
	//same as project to point that function is not used eleswhere so it has been eliminated
	double H[FEElement::MAX_NODES];
	double val = 0.0;
	//should be zero to proprly accumulate the values

	elem->shape_fnc(H, pos.x, pos.y, pos.z);
	for (size_t i = 0; i < elem->m_node.size(); i++)
	{
		val += gx[i] * H[i];
	}
	return val;
}

double FEAngio::FindECMDensity(const GridPoint& pt)
{
	assert(pt.nelem != -1 && pt.nelem != 0);
	FEMesh & mesh = m_fem->GetMesh();

	//TODO: replace with a check if it is in the same element before doing a search
	//the element may change between accesses
	/*
	double rez[3];
	FESolidElement* se = mesh.FindSolidElement(pt.r, rez);//TODO uses spatial coordinates
	assert(se->GetID_ang() == pt.nelem);//verify that migration is not happening
	*/
	FESolidElement * se = nullptr;
	if (pt.elemindex >= 0)
		se = dynamic_cast<FESolidElement*>(&pt.ndomain->ElementRef(pt.elemindex));
	else
		assert(false);
	double rez[3];
	rez[0] = pt.q.x; rez[1] = pt.q.y; rez[2] = pt.q.z;

	if (se)
	{
		//double * shapef = new double[se->Nodes()];
		double shapef[FEElement::MAX_NODES];
		se->shape_fnc(shapef, rez[0], rez[1], rez[2]);

		double coll_den = 0.0;
		for (int i = 0; i < se->Nodes(); i++)
		{
			int nn = se->m_node[i];
			nn = mesh.Node(nn).GetID();
			coll_den += m_fe_node_data[nn].m_ecm_den* shapef[i];
		}

		//delete[] shapef;
		return coll_den;
	}
	return 0.0;
}

//-----------------------------------------------------------------------------
void FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;
	FETimeInfo fti = fem.GetTime();
	m_time.t = fti.currentTime;
	m_time.dt = fti.timeIncrement;
	
	static bool start = false;
	if (!start)
	{
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->Update();
		}
		start = true;
	}

	if (nwhen == CB_UPDATE_TIME)
	{
		static int index = 0;
		// grab the time information
		
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->UpdateGDMs();
		}

		std::vector<std::pair<double,double>> times;
		index = FindGrowTimes(times, index);

		//new function to find the start time grow time and if this is the final iteration this timestep
		for (size_t i = 0; i < times.size(); i++)
		{
			FragmentBranching::Grow(times[i].first, times[i].second); //new grow method
		}
		

		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->AdjustMeshStiffness(m_pmat[i]->GetMaterial());
		}
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->UpdateSproutStressScaling();
		}
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->UpdateSprouts(1.0, m_pmat[i]->GetMatrixMaterial()->GetElasticMaterial());
		}
	}
	else if (nwhen == CB_MAJOR_ITERS)
	{
		// update the grid data
		UpdateECM();
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->Update();
		}

		++FE_state;
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			// Save the current vessel state
			fileout->save_vessel_state(*this);

			// save active tips
			fileout->save_active_tips(*this);

			// Print the status of angio3d to the user    
			fileout->printStatus(*this);
		}
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		if (!m_fem->GetGlobalConstant("no_io"))
		{
			Output();
		}
		//force any destructors to be called that need it
		delete fileout;
		fileout = nullptr;
			
	}
	//needed to copy data between material points to fix things for multiphasic materials
	else if(nwhen == CB_STEP_ACTIVE)
	{
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			//m_pmat[i]->ActiveFix();
		}
	}
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{	
	//write out the timeline of branchpoints
	fileout->save_timeline(*this);
	//fileout->save_winfiber(*this);
	fileout->save_final_vessel_csv(*this);
}

///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////
void FEAngio::adjust_mesh_stiffness()
{
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->AdjustMeshStiffness(m_pmat[i]->GetMaterial());
	}
}

//-----------------------------------------------------------------------------
void FEAngio::update_sprout_stress_scaling()
{
	//TODO: currently unused but needs moved or removed
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->scale = y0 + a / (1 + exp(-(m_time.t - x0) / b));
	}
		
	return;
}

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the fiber array
	int N = mesh.Nodes();
	fiber.assign(N, vec3d(0,0,0));

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local fiber orientation at integration points
			vector<double> fx(nint), fy(nint), fz(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mpoint->ExtractData<FEElasticMaterialPoint>();
				mat3d m = pt.m_Q;
				
				// grab the first column as the fiber orientation
				fx[n] = m[0][0];
				fy[n] = m[1][0];
				fz[n] = m[2][0];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln), gy(neln), gz(neln);
			el.project_to_nodes(&fx[0], &gx[0]);
			el.project_to_nodes(&fy[0], &gy[0]);
			el.project_to_nodes(&fz[0], &gz[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				vec3d f(gx[i], gy[i], gz[i]);
				fiber[ni] += f;
			}
		}
	}

	// normalize the fibers
	for (int i=0; i<N; ++i) fiber[i].unit();

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	density.resize(N, 0.0);
	anisotropy.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> den(nint);
			vector<double> anis(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
				den[n] = angioPt->m_D;
				anis[n] = angioPt->m_DA;
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			vector<double> fx(neln);
			el.project_to_nodes(&den[0], &gx[0]);
			el.project_to_nodes(&anis[0], &fx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				density[ni] = gx[i];
				if(fx[i] > anisotropy[ni])
					anisotropy[ni] = fx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	concentration.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> con(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FESolutesMaterialPoint& spt = *mpoint->ExtractData<FESolutesMaterialPoint>();
				con[n]=spt.m_c[vegfID];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			el.project_to_nodes(&con[0], &gx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				concentration[ni] = gx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}