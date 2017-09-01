#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include "Elem.h"
#include "Culture.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEFiberMaterialPoint.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FEBioMech/FESPRProjection.h"
#include <iostream>
#include "angio3d.h"


const double PI = 3.141592653589793;

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	ADD_PARAMETER2(m_cultureParams.sprout_s_mag, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "a");
	ADD_PARAMETER2(m_cultureParams.sprout_s_range, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "b");
	ADD_PARAMETER2(m_cultureParams.sprout_s_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "N");

	ADD_PARAMETER2(m_cultureParams.m_length_adjustment, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "length_adjustment");
	ADD_PARAMETER2(m_cultureParams.m_vessel_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "vessel_width");
	ADD_PARAMETER(m_cultureParams.growth_length_over_time, FE_PARAM_DOUBLE, "growth_length_over_time");

	ADD_PARAMETER(m_cultureParams.ecm_control, FE_PARAM_INT, "ecm_seeder");
	ADD_PARAMETER2(m_cultureParams.m_matrix_density, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "matrix_density");

	ADD_PARAMETER(m_cultureParams.m_symmetry_plane, FE_PARAM_VEC3D, "symmetryplane");
	//uncategorized variables are incomplete
	ADD_PARAMETER(m_cultureParams.m_composite_material, FE_PARAM_INT, "composite_material");
	ADD_PARAMETER2(m_cultureParams.m_sprout_force, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "sprout_force");

	ADD_PARAMETER2(m_cultureParams.active_tip_threshold, FE_PARAM_INT, FE_RANGE_GREATER_OR_EQUAL(0), "active_tip_threshold");
	ADD_PARAMETER2(m_cultureParams.stress_radius, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "stress_radius");
	
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticFiberMaterial(pfem), FEAngioMaterialBase()
{
	AddProperty(&common_properties, "common_properties");
	AddProperty(&matrix_material, "matrix");
}

FEAngioMaterial::~FEAngioMaterial()
{
	if (m_cult)
	{
		delete m_cult;
	}
	m_cult = nullptr;
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	// Create symmetry vectors
	sym_planes[0] = 0; sym_planes[1] = 0; sym_planes[2] = 0; sym_planes[3] = 0; sym_planes[4] = 0; sym_planes[5] = 0; sym_planes[6] = 0;

	sym_vects[0][0] = 1.; sym_vects[0][1] = 0.; sym_vects[0][2] = 0.;
	sym_vects[1][0] = 0.; sym_vects[1][1] = 1.; sym_vects[1][2] = 0.;
	sym_vects[2][0] = 0.; sym_vects[2][1] = 0.; sym_vects[2][2] = 1.;
	sym_vects[3][0] = 1.; sym_vects[3][1] = 1.; sym_vects[3][2] = 0.;
	sym_vects[4][0] = 1.; sym_vects[4][1] = 0.; sym_vects[4][2] = 1.;
	sym_vects[5][0] = 0.; sym_vects[5][1] = 1.; sym_vects[5][2] = 1.;
	sym_vects[6][0] = 1.; sym_vects[6][1] = 1.; sym_vects[6][2] = 1.;

	sym_on = false;

	if(!matrix_material->Init()) return false;

	if(!common_properties->vessel_material->Init()) return false;

	if (FEElasticMaterial::Init() == false) return false;

	//culture must be initialized here  so pangio is defined
	assert(m_pangio);
	m_cult = new Culture(*m_pangio, this, &m_cultureParams, common_properties->fbrancher);

	switch (m_cultureParams.ecm_control)
	{
	case 0:
		ecm_initializer = new ECMInitializerConstant();
		break;
	case 1:
		ecm_initializer = new ECMInitializerSpecified();
		break;
	case 2:
		ecm_initializer = new ECMInitializerNoOverwrite();
		break;
	default:
		assert(false);
	}

	// add the user sprouts
	std::vector<int> matls;
	
	FEMesh& mesh = GetFEModel()->GetMesh();
	
	FECoreBase * base = GetParent();
	if(base)
	{
		matls.emplace_back(base->GetID());
	}
	else
	{
		matls.emplace_back(this->GetID_ang());
	}

	mesh.DomainListFromMaterial(matls, domains);
	for (size_t i = 0; i < domains.size(); i++)
	{
		domainptrs.emplace_back(&mesh.Domain(domains[i]));
	}
	int co = 0;
	for (auto i = 0; i < mesh.Domains(); i++)
	{
		if (std::find(domainptrs.begin(), domainptrs.end(),&mesh.Domain(i)) != domainptrs.end())
		{
			meshOffsets[&mesh.Domain(i)] = co;
		}
		co += mesh.Domain(i).Elements();
	}
	for (unsigned int i=0; i<m_suser.size(); ++i)
	{
		if (domains.size())
			AddSprout(m_suser[i], vec3d(0,0,0), &mesh.Domain(domains[0]), matrix_material->GetElasticMaterial());
		//TODO: sprouts probably need distributed among the domains of the material
	}
	m_suser.clear();

	fiber_manager = new FiberManager(this);

	return true;
}
void FEAngioMaterial::FinalizeInit()
{
	FEMesh * mesh = m_pangio->GetMesh();
	// initialize material point data
	vec3d x[FEElement::MAX_NODES];

	for (int n = 0; n<mesh->Domains(); ++n)
	{
		FESolidDomain& dom = reinterpret_cast<FESolidDomain&>(mesh->Domain(n));
		FEMaterial* pm = dom.GetMaterial();
		FEAngioMaterial* pam = dynamic_cast<FEAngioMaterial*>(pm);

		if(!pam)
		{
			FEMultiphasic * mmat = dynamic_cast<FEMultiphasic*>(pm);
			if(mmat)
			{
				pam = dynamic_cast<FEAngioMaterial*>(mmat->GetSolid());
			}
			//consider adding the rest of the phasic materials
		}

		if (pam == this)
		{
			// loop over all elements
			int NE = dom.Elements();
			for (int i = 0; i<NE; ++i)
			{
				// get the next element
				FEElement& el = dom.Element(i);
				int neln = el.Nodes();

				// get the nodal coordinates
				for (int j = 0; j<neln; ++j) x[j] = mesh->Node(el.m_node[j]).m_rt;

				// loop over all integration points
				int nint = el.GaussPoints();
				for (int j = 0; j<nint; ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					FEAngioMaterialPoint* pt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
					if (pt)
					{
						vec3d r = el.Evaluate(x, j);

						// calculate the GridPoint data for this point.
						//TODO: check elastic material for integration point coordinates
						if (FindGridPoint(r, &dom, i, pt->m_pt) == false)
						{
							assert(false);
						}
					}
				}
			}
		}
		else
		{
			assert(false);
		}
	}
	for (int i = 0; i < domainptrs[0]->Nodes(); i++)
	{
		int nn = domainptrs[0]->NodeIndex(i);
		node_map[nn] = i;
	}
}
void FEAngioMaterial::SetupSurface()
{
	FEMesh * mesh = m_pangio->GetMesh();
	
	//setup the exterior_surface
	assert(domainptrs.size());
	exterior_surface = mesh->ElementBoundarySurface(domainptrs, true, false);
	normal_proj = new FENormalProjection(*exterior_surface);
	normal_proj->SetTolerance(0.001);
	normal_proj->Init();
	normal_proj->SetSearchRadius(1.0);

	//now add the exterior_surface element indices to the element data
	for (auto i = 0; i < exterior_surface->Elements(); i++)
	{
		FESurfaceElement & surfe = exterior_surface->Element(i);
		auto base_eindex = surfe.m_elem[0];

		m_pangio->m_fe_element_data[base_eindex + 1].surfacesIndices.emplace_back(i);
	}
}

void FEAngioMaterial::InitializeFibers()
{
	common_properties->InitializeFibers(fiber_manager);
}

void FEAngioMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	// get the material's coordinate system (if defined)
	FECoordSysMap* pmap = GetCoordinateSystemMap();
	//this allows the local coordinates to work correctly
	if (pmap)
	{
		FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
		FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
		FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);

		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		// compound the local map with the global material axes
		//mat3d Qlocal = pmap->LocalElementCoord(el, n);
		//pt.m_Q = pt.m_Q * Qlocal;

		vessel_elastic.m_Q = pt.m_Q;
		matrix_elastic.m_Q = pt.m_Q;

		FEElasticMaterial* vess_elastic = common_properties->vessel_material->GetElasticMaterial();
		FEElasticMaterial* mat_elastic = matrix_material->GetElasticMaterial();

		vess_elastic->SetLocalCoordinateSystem(el, n, *angioPt->vessPt);
		mat_elastic->SetLocalCoordinateSystem(el, n, *angioPt->matPt);
	}
	
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.name(), "sprout") == 0)
	{
		m_suser.emplace_back(m_s);
	}
}


//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::AngioStress(FEAngioMaterialPoint& angioPt)
{
	mat3ds s;
	s.zero();

	// get density scale factor
	double den_scale = m_cult->FindDensityScale(angioPt.m_pt);

	// loop over all sprout tips
	int NS = Sprouts();
	//TODO: fix the stress analysis
	
	// current position of integration point
	FESolidDomain * d = angioPt.m_pt.ndomain;

	vec3d y;
	assert(angioPt.m_pt.elemindex >= 0);
	y = CurrentPosition(&d->Element(angioPt.m_pt.elemindex), angioPt.m_pt.q.x, angioPt.m_pt.q.x, angioPt.m_pt.q.x);
		
	if (sym_on)
	{
//#pragma omp parallel for shared(s)
		for (int i = 0; i<NS; ++i)
		{
			SPROUT& sp = m_spr[i];

			// current position of sprout force
			vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);

			vec3d r = y - x;
			double l = r.unit();

			sp.sprout.unit();															// Normalize the sprout direction vector

			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			//TODO: some of this may be precalculated
			double p = den_scale*scale*m_cultureParams.sprout_s_mag*(pow(cos(theta / 2), m_cultureParams.sprout_s_width))*exp(-m_cultureParams.sprout_s_range*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

			//do we care about the sign of p?

			mat3ds si = dyad(r)*p;
														// If symmetry is turned on, apply symmetry
			MirrorSym(y, si, sp, den_scale);

//#pragma omp critical
			s += si;
		}
	}
	else
	{
		if (NS <= m_cultureParams.active_tip_threshold)
		{
//#pragma omp parallel for schedule(dynamic, 24)
			for (int i = 0; i<NS; ++i)
			{
				SPROUT& sp = m_spr[i];

				// current position of sprout force
				vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);

				vec3d r = y - x;
				double l = r.unit();

				sp.sprout.unit();															// Normalize the sprout direction vector

				double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

				//TODO: some of this may be precalculated
				double p = den_scale*scale*m_cultureParams.sprout_s_mag*(pow(cos(theta / 2), m_cultureParams.sprout_s_width))*exp(-m_cultureParams.sprout_s_range*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

				//double p = sprout_s_mag*exp(-sprout_s_range*l);

				mat3ds si = dyad(r)*p;

//#pragma omp critical
				s += si;
			}
		}
		else
		{
			std::vector<SPROUT>	temp;
			double local[3];
			local[0] = angioPt.m_pt.q.x;
			local[1] = angioPt.m_pt.q.y;
			local[2] = angioPt.m_pt.q.z;
			temp.emplace_back(vec3d(), &angioPt.m_pt.ndomain->Element(angioPt.m_pt.elemindex), local, dynamic_cast<FEAngioMaterialBase*>(this), matrix_material->GetElasticMaterial());
			std::pair<size_t, std::vector<SPROUT> *> dim = std::pair<size_t, std::vector<SPROUT> * >(0, &temp);
			std::vector<std::pair<size_t, std::vector<SPROUT> *>> nst;
			sprouts.within(dim, m_cultureParams.stress_radius * m_cultureParams.stress_radius, nst);
//#pragma omp parallel for schedule(dynamic, 24)
			for (int i = 0; i<nst.size(); ++i)
			{
				SPROUT& sp = m_spr[nst[i].first];

				// current position of sprout force
				vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);

				vec3d r = y - x;
				double l = r.unit();

				sp.sprout.unit();															// Normalize the sprout direction vector

				double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

				//TODO: some of this may be precalculated
				double p = den_scale*scale*m_cultureParams.sprout_s_mag*(pow(cos(theta / 2), m_cultureParams.sprout_s_width))*exp(-m_cultureParams.sprout_s_range*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

				//double p = sprout_s_mag*exp(-sprout_s_range*l);

				mat3ds si = dyad(r)*p;

//#pragma omp critical
				s += si;
			}
		}
	}
	return s;
}

void FEAngioMaterial::UpdateECM()
{
	ecm_initializer->updateECMdensity(this);
}

//this function accumulates the the anistropy and ecm_density, n_tag is incremented to be used to take the average
bool FEAngioMaterial::InitECMDensity(FEAngio * angio)
{
	ecm_initializer->seedECMDensity(this);
	return true;
}



mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	s.zero();
	//should always be true but we should check
	assert(angioPt);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

		vessel_elastic.m_rt = elastic_pt.m_rt;//spatial position
		vessel_elastic.m_r0 = elastic_pt.m_r0;//material position
		vessel_elastic.m_F = elastic_pt.m_F;//deformation gradient
		vessel_elastic.m_J = elastic_pt.m_J;//determinate
		
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;

		mat3ds activeStress = AngioStress(*angioPt) ;
		vessel_elastic.m_s = common_properties->vessel_material->Stress(*angioPt->vessPt);
		matrix_elastic.m_s = matrix_material->GetElasticMaterial()->Stress(*angioPt->matPt);

		s = activeStress + angioPt->vessel_weight*vessel_elastic.m_s + angioPt->matrix_weight*matrix_elastic.m_s;
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEAngioMaterial::Tangent(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	tens4ds s(0.0);
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		s = angioPt->vessel_weight*common_properties->vessel_material->Tangent(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->Tangent(*angioPt->matPt);
	}
	return s;
}

double FEAngioMaterial::StrainEnergyDensity(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
    
	// calculate strain energy density
	double sed = 0.0;
	if(angioPt)
	{
		FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
		vessel_elastic.m_rt = elastic_pt.m_rt;
		vessel_elastic.m_r0 = elastic_pt.m_r0;
		vessel_elastic.m_F = elastic_pt.m_F;
		vessel_elastic.m_J = elastic_pt.m_J;
		FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
		matrix_elastic.m_rt = elastic_pt.m_rt;
		matrix_elastic.m_r0 = elastic_pt.m_r0;
		matrix_elastic.m_F = elastic_pt.m_F;
		matrix_elastic.m_J = elastic_pt.m_J;
		sed = angioPt->vessel_weight*common_properties->vessel_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->matPt);
	}
	return sed;
}

void FEAngioMaterial::UpdateGDMs()
{
	common_properties->UpdateGDMs();
}


///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - ApplySym
//      Determine if symmetry is turned on, if so create the symmetry vectors
///////////////////////////////////////////////////////////////////////
void FEAngioMaterial::ApplySym()
{
	if (m_cultureParams.m_symmetry_plane.x != 0)															// Turn on x symmetry
		sym_planes[0] = 1;

	if (m_cultureParams.m_symmetry_plane.y != 0)															// Turn on y symmetry
		sym_planes[1] = 1;

	if (m_cultureParams.m_symmetry_plane.z != 0)															// Turn on z symmetry
		sym_planes[2] = 1;

	if ((m_cultureParams.m_symmetry_plane.x != 0) && (m_cultureParams.m_symmetry_plane.y != 0))												// Turn on x and y symmetry
		sym_planes[3] = 1;

	if ((m_cultureParams.m_symmetry_plane.x != 0) && (m_cultureParams.m_symmetry_plane.z != 0))												// Turn on x and z symmetry
		sym_planes[4] = 1;

	if ((m_cultureParams.m_symmetry_plane.y != 0) && (m_cultureParams.m_symmetry_plane.z != 0))												// Turn on y and z symmetry
		sym_planes[5] = 1;

	if ((m_cultureParams.m_symmetry_plane.x != 0) && (m_cultureParams.m_symmetry_plane.y != 0) && (m_cultureParams.m_symmetry_plane.z != 0))								// Turn on x y and z symmetry
		sym_planes[6] = 1;
	
	if (sym_planes[0] + sym_planes[1] + sym_planes[2] + sym_planes[3] + sym_planes[4] + sym_planes[5] + sym_planes[6] != 0)
		sym_on = true;														

	return;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEAngioMaterial::CreateMaterialPointData()
{
	return new FEAngioMaterialPoint(new FEFiberMaterialPoint(FEElasticMaterial::CreateMaterialPointData()), common_properties->vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData());
}