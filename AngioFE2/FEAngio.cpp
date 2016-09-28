///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMix/FESolute.h"
#include "FEBioMix/FEMultiphasic.h"
#include <FECore/log.h>
#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include "FECore/FESolidDomain.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "Elem.h"
#include "angio3d.h"
#include <ctime>
#include <algorithm>
#include <cfloat> 


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, vector<double>& anisotropy, FEMaterial* pmat);

bool CreateConcentrationMap(vector<double>& concentration, FEMaterial* pmat, int vegfID);

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : m_fem(fem)
{
	// Body force counter
	total_bdyf = 0;
	
	FE_state = 0;

	m_ntime = 1;

	// initialize time stepping parameters
	m_time.dt = 0.25;

	// Input random seed number
	m_irseed = 0;
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
}

//-----------------------------------------------------------------------------
FEModel& FEAngio::GetFEModel() const
{
	return m_fem;
}

FEMesh * FEAngio::GetMesh() const
{
	return &m_fem.GetMesh();
}

//-----------------------------------------------------------------------------
// find the angio material component
FEAngioMaterial* FEAngio::FindAngioMaterial(FEMaterial* pm)
{
	FEMaterial* pmat;
	if(strcmp(pm->GetTypeStr(), "angio")==0)
	{
		pmat = pm;
	}
	else
	{
		pmat = pm->FindComponentByType("angio");
	}
	if (pmat)
	{
		FEAngioMaterial* pma = dynamic_cast<FEAngioMaterial*>(pmat);
		return pma;
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	// Init all the FE stuff
	//must be done first initializes material
	if (InitFEM() == false) return false;

	// Seed the random number generator based on the sum of the seeds of the materials ie set the seed only once in any material
	unsigned int posseed = 0;
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		posseed += m_pmat[i]->GetSeed();
	}
	if (!posseed)
	{
		posseed = static_cast<unsigned int>(time(0));
	}
	m_irseed = posseed;
	srand(m_irseed);

	// Print out the seed number for the random generator
	fileout.printrandseed(m_irseed);



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


	//setup the exterior_surface
	SetupSurface();

	// assign ECM densities to grid nodes
	//TODO: degree of anisotropy values will only get initialized if density is set to 0
	if (InitECMDensity() == false) return false;

	// assign collagen fibers to grid nodes
	if (InitCollagenFibers() == false) return false;

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

void FEAngio::SetupSurface()
{
	for (int i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->SetupSurface();
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
	for (int i = 0; i < m_fem.Materials(); i++)
	{
		FEAngioMaterial * cmat = dynamic_cast<FEAngioMaterial*>(m_fem.GetMaterial(i));
		if (cmat)
		{
			m_pmat.push_back(cmat);
			m_pmat_ids.push_back(cmat->GetID());
			//TODO: check that material parameters are set here
			cmat->ApplySym();
			cmat->SetFEAngio(this);
			cmat->UpdateSproutStressScaling();
		}
	}
	assert(m_pmat.size());

	felog.printf("%d Angio materials found. Stress approach will be used.", m_pmat.size());

	// register the callback
	m_fem.AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED, this);

	// Do the model initialization
	if (m_fem.Init() == false) return false;

	for (int i = 0; i < m_fem.Materials(); i++)
	{
		FEAngioMaterial * cmat = dynamic_cast<FEAngioMaterial*>(m_fem.GetMaterial(i));
		if (cmat)
		{
			cmat->SetBoundaryCondition();
		}
	}
	return true;
}
void FEAngio::FinalizeFEM()
{
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		//TODO: remove this constant or make it a user parameter
		m_pmat[i]->CreateSprouts(0.5);
		m_pmat[i]->AdjustMeshStiffness();
	}

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::FILE_ONLY);

	// --- Output initial state of model ---

	// Output initial microvessel state
	fileout.save_vessel_state(*this);

	// save active tips
	fileout.save_active_tips(*this);

	// Output time information
	fileout.save_time(*this);

	// Output initial collagen fiber orientation
	fileout.writeCollFib(*this, true);
}

//-----------------------------------------------------------------------------
// Initialize the nodal ECM values
bool FEAngio::InitECMDensity()
{
	ForEachNode([this](FENode & node)
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
	ForEachNode([this](FENode & node)
	{
		//nneds to be run only once per node
		if (m_fe_node_data[node.GetID()].m_ntag)
		{
			m_fe_node_data[node.GetID()].m_ecm_den0 /= static_cast<double>(m_fe_node_data[node.GetID()].m_ntag);
			m_fe_node_data[node.GetID()].m_ecm_den = m_fe_node_data[node.GetID()].m_ecm_den0;
			//m_fe_node_data[node.GetID()].m_da /= static_cast<double>(m_fe_node_data[node.GetID()].m_ntag);
			m_fe_node_data[node.GetID()].m_ntag = 0;
		}
	});
	return true;
}
//-----------------------------------------------------------------------------
// update the extracellular matrix nodal values
void FEAngio::UpdateECM()
{
	//straight translation of code from the grid probably should be optimized later

	// reset nodal data
	FEMesh & mesh = m_fem.GetMesh();

	ForEachNode([this](FENode & node)
	{
		m_fe_node_data[node.GetID()].m_collfib = vec3d(0, 0, 0);
		m_fe_node_data[node.GetID()].m_ecm_den = 0.0;
		//REFACTOR: why reset not just overwrite
	});

	//this portion will be harder
	// For each element within the grid...
	ForEachElement([this, &mesh](FESolidElement & elem, FESolidDomain & d)
	{
		//these will hold the natural coordinates once the project to nodes is complete 
		double nr[FEElement::MAX_NODES];
		double ns[FEElement::MAX_NODES];
		double nt[FEElement::MAX_NODES];
		//these hold the natural coordinates of the integration points (r,s,t)
		double gr[FEElement::MAX_NODES];
		double gs[FEElement::MAX_NODES];
		double gt[FEElement::MAX_NODES];
		//TODO: if needed get FEBIO to expose the vectors that contain these to avoid this copy
		for (int i = 0; i < elem.Nodes(); i++)
		{
			gr[i] = elem.gr(i);
			gs[i] = elem.gs(i);
			gt[i] = elem.gt(i);
		}

		elem.project_to_nodes(gr, nr);
		elem.project_to_nodes(gs, ns);
		elem.project_to_nodes(gt, nt);

		// For each node in the element...
		for (int j = 0; j<elem.Nodes(); ++j)
		{
			// get the node
			int nnum = elem.m_node[j];
			nnum = mesh.Node(nnum).GetID();
			// get the ecm density and collagen fiber
			double ecm_den = m_fe_node_data[nnum].m_ecm_den0;
			vec3d coll_fib = m_fe_node_data[nnum].m_collfib0;

			/*
			//clamp n* to [1,-1]
			nr[j] = min(max(nr[j], -1), 1);
			ns[j] = min(max(ns[j], -1), 1);
			nt[j] = min(max(nt[j], -1), 1);
			*/

			//round to nearest integer
			nr[j] = round(nr[j]);
			ns[j] = round(ns[j]);
			nt[j] = round(nt[j]);

			// Calculate the deformation gradient tensor and jacobian at the node
			mat3d F;
			double Jacob = d.defgrad(elem, F, nr[j], ns[j], nt[j]);
			
			//make sure the function is differentiable and preserves orientation
			assert(Jacob > 0.0);

			// Update the collagen fiber orientation vector into the current configuration using F		
			coll_fib = F*coll_fib;
			coll_fib.unit();

			// Update matrix density using the Jacobian
			ecm_den = ecm_den / Jacob;

			// accumulate fiber directions and densities
			m_fe_node_data[nnum].m_collfib += coll_fib;
			m_fe_node_data[nnum].m_ecm_den += ecm_den;


			// increment counter
			m_fe_node_data[nnum].m_ntag++;
		}
	  });

	// normalize fiber vector and average ecm density
	ForEachNode([this](FENode & node)
	{
		//nneds to be run only once per node
		if (m_fe_node_data[node.GetID()].m_ntag)
		{
			m_fe_node_data[node.GetID()].m_ecm_den /= (double)m_fe_node_data[node.GetID()].m_ntag;
			m_fe_node_data[node.GetID()].m_collfib.unit();
			m_fe_node_data[node.GetID()].m_ntag = 0;
		}
	});
}

//-----------------------------------------------------------------------------
// Find VEGF solute
int FEAngio::FindVEGF()
{
	FEModel& fem = GetFEModel();
	int N = fem.GlobalDataItems();
	for (int i=0; i<N; ++i)
	{
		FESoluteData* psd = dynamic_cast<FESoluteData*>(fem.GetGlobalData(i));
		if (psd && (strcmp(psd->m_szname,"VEGF"))) return psd->m_nID;
	}
	return 0;
}
//-----------------------------------------------------------------------------
// Initialize the nodal concentration values
bool FEAngio::InitSoluteConcentration()
{
	FEMesh & mesh = m_fem.GetMesh();
	int vegfID = FindVEGF();

	if(vegfID==0)
		return false;

	int NN = mesh.Nodes();
	vector<double> concentration(NN, 0.0);

	// get the material
	FEMaterial* pm = m_fem.GetMaterial(0);
	if (CreateConcentrationMap(concentration, pm, vegfID) == false) return false;

	// assign ECM density
	for (int i = 0; i < NN; ++i)								
	{
		FEAngioNodeData& node = m_fe_node_data[i +1];

		node.vegf_conc = concentration[i];
	}

	return true;
}

//TODO: consider making the ForEach* const as long as this works on both compilers
void FEAngio::ForEachNode(std::function<void(FENode &)> f, std::vector<int> & matls)
{
	//TODO: the last element to access a node wins on overwting the data associated with that node
	//this behavior matches the previous behavior of the plugin but probably should be fixed sometime
	FEMesh & mesh = m_fem.GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	for (size_t i = 0; i < dl.size(); i++)
	{
		FEDomain & d = mesh.Domain(dl[i]);
		for (int j = 0; j < d.Elements(); j++)
		{
			FEElement & e = d.ElementRef(j);
			for (int k = 0; k < e.Nodes(); k++)
			{
				f(mesh.Node(e.m_node[k]));//this iterates over the local nodes
			}
		}
	}
}
void FEAngio::ForEachNode(std::function<void(FENode &)> f)
{
	ForEachNode(f, m_pmat_ids);
}
void FEAngio::ForEachElement(std::function<void(FESolidElement&, FESolidDomain&)> f, std::vector<int> & matls)
{
	FEMesh & mesh = m_fem.GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	for (size_t i = 0; i < dl.size(); i++)
	{
		FESolidDomain & d = reinterpret_cast<FESolidDomain&>(mesh.Domain(i));
		for (int j = 0; j < d.Elements(); j++)
		{
			FESolidElement & e = reinterpret_cast<FESolidElement&>(d.ElementRef(j));
			f(e, d);
		}
	}
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
	FEMesh & mesh = m_fem.GetMesh();
	std::vector<int> dl;
	mesh.DomainListFromMaterial(matls, dl);
	for (size_t i = 0; i < dl.size(); i++)
	{
		FESolidDomain & d = reinterpret_cast<FESolidDomain&>(mesh.Domain(i));
		f(d);
	}
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
	size_t j;
	// get the element nodal coordinates
	auto neln = se->Nodes();
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
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(1, rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//s
	cpos = vec3d(rs[0], 1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], 1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[0], -1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], -1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//t
	cpos = vec3d(rs[0], rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[0], rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//do all of the negations
	//first
	//r
	cpos = vec3d(1, -rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(1, -rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, -rs[0], rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, -rs[1], rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//s
	cpos = vec3d(-rs[0], 1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], 1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[0], -1, rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], -1, rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//t
	cpos = vec3d(-rs[0], rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[0], rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//second
	//r
	cpos = vec3d(1, rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(1, rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//s
	cpos = vec3d(rs[0], 1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], 1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[0], -1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], -1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//t
	cpos = vec3d(rs[0], -rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], -rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[0], -rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(rs[1], -rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);


	//both
	//r
	cpos = vec3d(1, -rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(1, -rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, -rs[0], -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-1, -rs[1], -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//s
	cpos = vec3d(-rs[0], 1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], 1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[0], -1, -rs[1]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], -1, -rs[0]);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	//t
	cpos = vec3d(-rs[0], -rs[1], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], -rs[0], 1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[0], -rs[1], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);

	cpos = vec3d(-rs[1], -rs[0], -1);
	cur = LocalToGlobal(elem, cpos);
	positions_global.push_back(cur);
	positions_local.push_back(cpos);


	double best_dist = DBL_MAX;
	int best_index = -1;
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

GridPoint FEAngio::FindGridPoint(FEDomain * domain, int nelem, vec3d& q) const
{
	assert(domain != nullptr && nelem >= 0);
	GridPoint pt;
	pt.q = q;
	FEMesh & mesh = m_fem.GetMesh();
	FESolidElement * se;
	//TODO: refactor if problems with multiple domains
	if ((se = dynamic_cast<FESolidElement*>(&domain->ElementRef(nelem))))
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
	FEMesh & mesh = m_fem.GetMesh();
	vec3d r(0, 0, 0);
	FEDomain * d = pt.ndomain;
	FESolidElement * se;
	if ((se = dynamic_cast<FESolidElement*>(&d->ElementRef(pt.elemindex))))
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

std::vector<double> FEAngio::createVectorOfMaterialParameters(FEDomain & d, FEElement * elem,
	double FEAngioNodeData::*materialparam, double r, double s, double t)
{
	FEMesh * mesh = GetMesh();
	std::vector<double> gx(elem->m_node.size());
	for (size_t i = 0; i < elem->m_node.size(); i++)
	{
		gx[i] = this->m_fe_node_data[mesh->Node(elem->m_node[i]).GetID()].*materialparam;
	}
	return gx;
}
double FEAngio::genericProjectToPoint(FEDomain & d, FEElement * elem,
	double FEAngioNodeData::*materialparam, double r, double s, double t)
{
	std::vector<double> gx = createVectorOfMaterialParameters(d, elem, materialparam, r, s, t);
	//same as project to point that function is not used eleswhere so it has been eliminated
	double H[FEElement::MAX_NODES];
	double val = 0.0;
	//should be zero to proprly accumulate the values

	//check that 0 if not solid element is okay
	FESolidElement * se = dynamic_cast<FESolidElement*>(elem);
	if (se)
	{
		se->shape_fnc(H, r, s, t);
		for (size_t i = 0; i < elem->m_node.size(); i++)
		{
			val += gx[i] * H[i];
		}
	}
	return val;
}

double FEAngio::FindECMDensity(const GridPoint& pt)
{
	assert(pt.nelem != -1 && pt.nelem != 0);
	FEMesh & mesh = m_fem.GetMesh();

	//TODO: replace with a check if it is in the same element before doing a search
	//the element may change between accesses
	/*
	double rez[3];
	FESolidElement* se = mesh.FindSolidElement(pt.r, rez);//TODO uses spatial coordinates
	assert(se->GetID() == pt.nelem);//verify that migration is not happening
	*/
	FESolidElement * se;
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

vec3d FEAngio::CollagenDirection(GridPoint& pt)
{
	assert(pt.elemindex >= 0);
	assert(pt.ndomain != nullptr);
	FEMesh & mesh = m_fem.GetMesh();
	// get the element
	FEElement * elem = &pt.ndomain->ElementRef(pt.elemindex);

	FESolidElement * selem = dynamic_cast<FESolidElement *>(elem);
	
	//TODO: may be refactored to remove shape function and dependencies 
	// Obtain shape function weights
	double shapeF[FEElement::MAX_NODES];
	if (selem)
		selem->shape_fnc(shapeF, pt.q.x, pt.q.y, pt.q.z);//check these numbers

	// Determine component of new vessel direction due to nodal collagen fiber orientation
	vec3d coll_angle(0, 0, 0);
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		coll_angle += m_fe_node_data[mesh.Node(n).GetID()].m_collfib*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}

//-----------------------------------------------------------------------------
// Initialize nodal collagen fiber directions
bool FEAngio::InitCollagenFibers()
{
	bool rv = true;
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		rv &= m_pmat[i]->InitCollagenFibers();
	}
	return rv;
}

//-----------------------------------------------------------------------------
void FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;
	m_time.t = fem.m_ftime;

	if (nwhen == CB_UPDATE_TIME)
	{
		// grab the time information
		
		m_time.dt = fem.GetCurrentStep()->m_dt;

		// do a growth step
		for (size_t i = 0; i < m_pmat.size(); i++)
		{
			m_pmat[i]->Grow(m_time);
			m_pmat[i]->AdjustMeshStiffness();
			m_pmat[i]->UpdateSproutStressScaling();
			m_pmat[i]->UpdateSprouts(1.0);
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

		// Save the current vessel state
		fileout.save_vessel_state(*this);

		// save active tips
		fileout.save_active_tips(*this);

		// Output time information	
		fileout.save_time(*this);
		
		// Print the status of angio3d to the user    
		fileout.printStatus(*this);
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		Output();
	}
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{
	// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	fileout.output_params(*this);
	
	// Output data file
	fileout.dataout(*this);
		
	// Output final collagen fiber orientation
	fileout.writeCollFib(*this, false);

	// Output final matrix density
	fileout.writeECMDen(*this);
}

///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////
// TODO: vessel lengths are always positive now, so we need to fix the logic here.
void FEAngio::adjust_mesh_stiffness()
{
	for (size_t i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->AdjustMeshStiffness();
	}
}

//-----------------------------------------------------------------------------
void FEAngio::update_sprout_stress_scaling()
{
	//TODO: currently unused but needs moved or removed
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	for (int i = 0; i < m_pmat.size(); i++)
	{
		m_pmat[i]->scale = y0 + a / (1 + exp(-(m_time.t - x0) / b));
	}
		
	return;
}

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
	// get the material's coordinate system
	FECoordSysMap* pmap = pmat->GetCoordinateSystemMap();

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