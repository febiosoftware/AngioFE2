#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include "Elem.h"
#include "Culture.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include <iostream>
#include "angio3d.h"


const double PI = 3.141592653589793;

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEAngioMaterialPoint, FEMaterialPoint)
	ADD_PARAMETER(m_D, FE_PARAM_DOUBLE, "dens");
	ADD_PARAMETER(m_DA, FE_PARAM_DOUBLE, "anisotropy");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterialPoint::FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt) : FEMaterialPoint(pt)
{
	m_D = 0.0;
	m_DA = 1.0;
	vessPt = vesselPt;
	matPt = matrixPt;
	vessPt->SetPrev(this);
	matPt->SetPrev(this);
	m_D = 0.0;
}


//-----------------------------------------------------------------------------
//! The init function is used to intialize data
void FEAngioMaterialPoint::Init()
{
	FEMaterialPoint::Init();
	vessPt->Init();
	matPt->Init();
}

//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
FEMaterialPoint* FEAngioMaterialPoint::Copy()
{
	FEAngioMaterialPoint* pt = new FEAngioMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
void FEAngioMaterialPoint::Serialize(DumpStream& dmp)
{
	if (dmp.IsSaving())
	{
		dmp << m_D;
	}
	else
	{
		dmp >> m_D;
	}
	FEMaterialPoint::Serialize(dmp);
}

FEAngioMaterialPoint* FEAngioMaterialPoint::FindAngioMaterialPoint(FEMaterialPoint* mp)
{
	FEAngioMaterialPoint* angioPt  = dynamic_cast<FEAngioMaterialPoint*>(mp);
	if(angioPt)
		return angioPt;

	FEMaterialPoint* pt = mp;
	while(pt)
	{
		angioPt = dynamic_cast<FEAngioMaterialPoint*>(pt);
		if(angioPt)
			return angioPt;

		FEElasticMixtureMaterialPoint* mixtureP = dynamic_cast<FEElasticMixtureMaterialPoint*>(pt);
		if(mixtureP)
		{
			vector<FEMaterialPoint*> mixtureVector = mixtureP->m_mp;
			for(unsigned int i=0; i<mixtureVector.size(); i++)
			{
				//TODO: is the recursion needed or not(is this search too deep?)
				angioPt = FindAngioMaterialPoint(mixtureVector[i]);
				if(angioPt)
				{
					return angioPt;
				}
			}
		}

		pt = pt->Next();
	}

	return nullptr;
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	ADD_PARAMETER2(m_cultureParams.sprout_s_mag, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "a");
	ADD_PARAMETER2(m_cultureParams.sprout_s_range, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "b");
	ADD_PARAMETER2(m_cultureParams.sprout_s_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "N");

	ADD_PARAMETER2(m_cultureParams.m_length_adjustment, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "length_adjustment");
	ADD_PARAMETER2(m_cultureParams.m_vessel_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "vessel_width");
	ADD_PARAMETER(m_cultureParams.growth_length_over_time, FE_PARAM_DOUBLE, "growth_length_over_time");

	ADD_PARAMETER2(m_cultureParams.m_matrix_condition, FE_PARAM_INT, FE_RANGE_CLOSED(0,4), "matrix_condition");
	ADD_PARAMETER(m_cultureParams.ecm_control, FE_PARAM_INT, "ecm_seeder");
	ADD_PARAMETER2(m_cultureParams.m_matrix_density, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "matrix_density");

	ADD_PARAMETER(m_cultureParams.m_symmetry_plane, FE_PARAM_VEC3D, "symmetryplane");
	//uncategorized variables are incomplete
	ADD_PARAMETER(m_cultureParams.m_composite_material, FE_PARAM_INT, "composite_material");
	ADD_PARAMETER2(m_cultureParams.m_sprout_force, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "sprout_force");

	ADD_PARAMETER2(m_cultureParams.active_tip_threshold, FE_PARAM_INT, FE_RANGE_GREATER_OR_EQUAL(0), "active_tip_threshold");
	ADD_PARAMETER2(m_cultureParams.stress_radius, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "stress_radius");
	
END_PARAMETER_LIST();

FEAngioMaterial::SPROUT::SPROUT(vec3d dir, FEElement * el, double * local, FEAngioMaterial * m) : sprout(dir), pel(el), mat(m)
{
	r[0] = local[0];
	r[1] = local[1];
	r[2] = local[2];
}

std::vector<double> units3d(3, 1.0);

std::vector<double> access_sprout(std::pair<size_t, std::vector<FEAngioMaterial::SPROUT> *> p)
{
	std::vector<double> rv;
	FEAngioMaterial::SPROUT spr = (*p.second)[p.first];
	vec3d cpos = spr.mat->CurrentPosition(spr.pel, spr.r[0], spr.r[1], spr.r[2]);
	rv.push_back(cpos.x);
	rv.push_back(cpos.y);
	rv.push_back(cpos.z);
	return rv;
}

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem), sprouts(access_sprout, ndim_distance, ndim_distance_to_plane, units3d)
{
	scale = 1.0;
	
	m_pangio = nullptr;

	m_cult = nullptr;

	AddProperty(&vessel_material, "vessel");
	AddProperty(&matrix_material , "matrix");
	AddProperty(&fbrancher, "brancher");
	AddProperty(&fseeder, "fragment_seeder");
	AddProperty(&bc, "boundary_condition");
	AddProperty(&gdms, "grow_direction_modifiers");

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

	if(!vessel_material->Init()) return false;

	if (FEElasticMaterial::Init() == false) return false;


	//culture must be initialized here  so pangio is defined
	assert(m_pangio);
	m_cult = new Culture(*m_pangio, this, &m_cultureParams, fbrancher);

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
	matls.push_back(this->GetID());
	FEMesh& mesh = GetFEModel()->GetMesh();
	mesh.DomainListFromMaterial(matls, domains);
	for (size_t i = 0; i < domains.size(); i++)
	{
		domainptrs.push_back(&mesh.Domain(domains[i]));
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
			AddSprout(m_suser[i], vec3d(0,0,0), &mesh.Domain(domains[0]));
		//TODO: sprouts probably need distributed among the domains of the material
	}
	m_suser.clear();

	

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
		FEAngioMaterial* pam;
		if (strcmp(pm->GetTypeStr(), "angio") == 0)
		{
			pam = dynamic_cast<FEAngioMaterial*>(pm);
		}
		else
		{
			assert(false);
			//pam = dynamic_cast<FEAngioMaterial*>(pm->FindComponentByType("angio"));
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

		m_pangio->m_fe_element_data[base_eindex + 1].surfacesIndices.push_back(i);
	}
}

vec3d FEAngioMaterial::CollagenDirection(GridPoint& pt)
{
	assert(pt.elemindex >= 0);
	assert(pt.ndomain != nullptr);
	FEMesh & mesh = m_pangio->m_fem.GetMesh();
	// get the element
	FEElement * elem = &pt.ndomain->ElementRef(pt.elemindex);

	FESolidElement * selem = dynamic_cast<FESolidElement *>(elem);

	//TODO: may be refactored to remove shape function and dependencies 
	// Obtain shape function weights
	double shapeF[FEElement::MAX_NODES];
	assert(selem);
	selem->shape_fnc(shapeF, pt.q.x, pt.q.y, pt.q.z);//check these numbers

	// Determine component of new vessel direction due to nodal collagen fiber orientation
	vec3d coll_angle(0, 0, 0);
	for (int i = 0; i<8; ++i)
	{
		int n = elem->m_node[i];
		coll_angle += m_pangio->m_fe_node_data[mesh.Node(n).GetID()].m_collfib*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

	return coll_angle;
}

void FEAngioMaterial::AdjustMeshStiffness()
{
	FEMesh & mesh = m_pangio->m_fem.GetMesh();
	if (m_cultureParams.m_composite_material == 0)													// If a composite consitutive model isn't being used, exit
		return;

	int elem_num = 0;													// Element number
	vec3d vess_vect;													// Vessel vector
	std::vector<int> matls;
	matls.push_back(GetID());
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1 / static_cast<double>(Nsub);									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

	//Zero the element items needed
	m_pangio->ForEachElement([this](FESolidElement & se, FESolidDomain & d)
	{
		int elemnum = se.GetID();
		m_pangio->m_fe_element_data[elemnum].alpha = 0.0;
		m_pangio->m_fe_element_data[elemnum].fiber_orient = vec3d(0, 0, 0);
	}, matls);

	const SegmentList& seg_list = m_cult->GetSegmentList();
	for (ConstSegIter frag_it = seg_list.begin(); frag_it != seg_list.end(); ++frag_it)		// For each segment...
	{
		Segment subunit;												// Segment subdivision placeholder

		const Segment& seg = (*frag_it);												// Obtain the segment

		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			assert(seg.length() > 0.);
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				subunit.tip(0).pt = seg.tip(0).pt;
			}

			// Calculate the subdivision
			vec3d v = seg.uvect();
			subunit.tip(1).pt.r = subunit.tip(0).pos() + v*(sub_scale*seg.length());

			subunit.Update();										// Find the length of the subdivision
			double intp[3];
			FEElement * volElement;
			mid = (subunit.tip(1).pos() + subunit.tip(0).pos())*0.5;
			if ((subunit.tip(0).pt.elemindex != -1) && 
				m_pangio->IsInsideHex8(&static_cast<FESolidElement&>(subunit.tip(0).pt.ndomain->ElementRef(subunit.tip(0).pt.elemindex)), mid, &mesh, intp))
			{
				elem_num = subunit.tip(0).pt.nelem;
				volElement = &(subunit.tip(0).pt.ndomain->ElementRef(subunit.tip(0).pt.elemindex));
			}
			else if ((subunit.tip(1).pt.elemindex != -1) &&
				m_pangio->IsInsideHex8(&static_cast<FESolidElement&>(subunit.tip(1).pt.ndomain->ElementRef(subunit.tip(1).pt.elemindex)),mid, &mesh, intp))
			{
				elem_num = subunit.tip(1).pt.nelem;
				volElement = &(subunit.tip(1).pt.ndomain->ElementRef(subunit.tip(1).pt.elemindex));
			}
			else
			{
				volElement = mesh.FindSolidElement(mid, intp);
				//TODO: be able to remove this check
				if (!volElement)
					continue;
				elem_num = volElement->GetID(); // Find the element that the midpoint is within
			}

			// Calculate the orientation of the subdivision
			vess_vect = subunit.tip(1).pos() - subunit.tip(0).pos();

			assert(vess_vect.norm() != 0.);									// Normalize to find the unit vector
			vess_vect = vess_vect / vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				FEAngioElementData& el = m_pangio->m_fe_element_data[elem_num];

				elem_volume = mesh.ElementVolume(*volElement);// Calculate the volume of the element

				subunit_volume = PI*(m_cultureParams.m_vessel_width / 2.)*(m_cultureParams.m_vessel_width / 2.)*fabs(subunit.length());		// Find the volume of the subdivision
				volume_fraction = subunit_volume / elem_volume;				// Calculate the volume fraction

				el.alpha = el.alpha + volume_fraction;	// Add the volume fraction for each subdivision to alpha

				// Calculate the vessel orientation vector 
				if ((el.fiber_orient.x == 0) && (el.fiber_orient.y == 0) && (el.fiber_orient.z == 0)){	// If the vessel orientation vector hasn't been assigned yet...
					el.fiber_orient = vess_vect;			// Set the vessel orientation vector					
				}
				else{														// If it has been...	
					el.fiber_orient.x = (el.fiber_orient.x + vess_vect.x) / 2;	// Average together the vessel orientation vector
					el.fiber_orient.y = (el.fiber_orient.y + vess_vect.y) / 2;
					el.fiber_orient.z = (el.fiber_orient.z + vess_vect.z) / 2;
				}
			}

			// Set the origin of the next subdivision to the end of the current one
			subunit.tip(0).pt.r = subunit.tip(1).pos();
		}
	}

	// Volume fraction for the composite material model




	m_pangio->ForEachElement([this, &mesh](FESolidElement & e, FESolidDomain & d)
	{
		assert(std::find(domainptrs.begin(), domainptrs.end(), &d) != domainptrs.end());
		vec3d e1; vec3d e2; vec3d e3;						// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
		double alpha = 0.;									// Obtain the element from the domain
		int nint = e.GaussPoints();										// Obtain the number of gauss points
		int num_elem = e.GetID();

		FEAngioElementData& eg = m_pangio->m_fe_element_data[num_elem];
		alpha = eg.alpha;											// Obtain alpha from the grid element
		for (int n = 0; n< e.Nodes(); n++)
		{
			int id = e.m_node[n];
			id = mesh.Node(id).GetID();
			m_pangio->m_fe_node_data[id].alpha = alpha;
		}

		// Set e1 to the vessel orientation vector
		e1 = eg.fiber_orient;

		if ((e1.x == 0) && (e1.y == 0) && (e1.z == 0)){						// If there is not vessels in the element, set the material basis to the global coordinate basis
			e1 = vec3d(1, 0, 0);
			e2 = vec3d(0, 1, 0);
			e3 = vec3d(0, 0, 1);
		}
		else{																// Else, set the other two directions to be orthogonal to the vessel orientation
			e2.y = 1;
			e2 = e1^e2;
			e3 = e1^e2;
		}

		for (int n = 0; n < nint; ++n)										// For each gauss point...
		{
			FEMaterialPoint& mp = *(e.GetMaterialPoint(n));
			FEAngioMaterialPoint* pt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp); // get the mixture material point
			pt->vessel_weight = alpha;
			pt->matrix_weight = 1.0 - alpha;

			if (m_cultureParams.m_composite_material == 2){													// If the transversely isotropic material is being used...
				FEElasticMaterialPoint& pt2 = *mp.ExtractData<FEElasticMaterialPoint>();
				pt2.m_Q[0][0] = e1.x;												// Set the first column of Q to e1
				pt2.m_Q[1][0] = e1.y;
				pt2.m_Q[2][0] = e1.z;
				pt2.m_Q[0][1] = e2.x;												// Set the second column of Q to e2
				pt2.m_Q[1][1] = e2.y;
				pt2.m_Q[2][1] = e2.z;
				pt2.m_Q[0][2] = e3.x;												// Set the third column of Q to e3
				pt2.m_Q[1][2] = e3.y;
				pt2.m_Q[2][2] = e3.z;
			}
		}

		num_elem++;


	}, matls);

	return;
}


bool FEAngioMaterial::FindGridPoint(const vec3d & r, GridPoint & p) const
{
	FEMesh * mesh = m_pangio->GetMesh();
	double natc[3];
	FESolidElement * se = mesh->FindSolidElement(r, natc);
	if (se && (std::find(domainptrs.begin(), domainptrs.end(), se->GetDomain()) != domainptrs.end()))
	{
		p.r = r;
		p.q.x = natc[0];
		p.q.y = natc[1];
		p.q.z = natc[2];
		p.ndomain = dynamic_cast<FESolidDomain*>(se->GetDomain());
		p.nelem = se->GetID();
		//TODO: hack
		p.elemindex = se->GetID() - 1 - meshOffsets.find(p.ndomain)->second;
		
		vec3d pq = m_pangio->Position(p);
		vec3d pw = p.r;
		assert((m_pangio->Position(p) - p.r).norm() < 1.0);
		return true;
	}
	p.ndomain = nullptr;
	p.nelem = -1;
	return false;
}

bool FEAngioMaterial::FindGridPoint(const vec3d & r, FESolidDomain * domain, int elemindex, GridPoint & p) const
{
	FEMesh * mesh = m_pangio->GetMesh();
	double natc[3];
	FESolidElement & se = reinterpret_cast<FESolidElement&>(domain->ElementRef(elemindex));
	if (m_pangio->IsInsideHex8(&se, r, mesh, natc))
	{
		p.r = r;
		p.q.x = natc[0];
		p.q.y = natc[1];
		p.q.z = natc[2];
		p.nelem = se.GetID();
		p.elemindex = elemindex;
		p.ndomain = domain;
		vec3d pq = m_pangio->Position(p);
		assert((m_pangio->Position(p) - p.r).norm() < 1.0);
		assert(se.GetDomain() == domain);
		return true;
	}
	return false;
}

bool FEAngioMaterial::InitCollagenFibers()
{
	std::vector<int> matls;
	matls.push_back(GetID());
	//remove later
	FEMesh * mesh = m_pangio->GetMesh();
	std::vector<vec3d> colfibs;
	colfibs.reserve(mesh->Nodes());
	//modes 1,3 are multimaterila safe
	switch (m_cultureParams.m_matrix_condition)
	{
	case 0: // random orientation
		//TODO: ask jeff if this should be removed?
		//this mode is the old way of doing things and should only be used for verifying that the test suite still passes
		for (int i = 0; i < mesh->Nodes(); i++)
		{
			vec3d v = vrand();
			if (m_cultureParams.m_bzfibflat) v.z *= 0.25;

			// normalize the vector
			v.unit();
			colfibs.push_back(v);
		}

		for (int i = 0; i < mesh->Nodes();i++)
		{
			// assign the node
			FENode & node = mesh->Node(i);
		
			vec3d v = colfibs[i];
			m_pangio->m_fe_node_data[node.GetID()].m_collfib0 = v;
			m_pangio->m_fe_node_data[node.GetID()].m_collfib = v;
		}//, matls);
		
		break;
	case 1: //multimaterial safe method

		m_pangio->ForEachNode([&](FENode & node)
		{
			vec3d v = m_pangio->uniformRandomDirection();
			v.unit();
			m_pangio->m_fe_node_data[node.GetID()].m_collfib0 = v;
			m_pangio->m_fe_node_data[node.GetID()].m_collfib = v;
		}, matls);
		break;
	case 3:	// from element's local coordinates
		m_pangio->ForEachElement([&](FESolidElement & se, FESolidDomain & sd)
		{
			size_t neln = se.Nodes();
			size_t nint = se.GaussPoints();

			// local fiber orientation at integration points
			double fx[FEElement::MAX_INTPOINTS], fy[FEElement::MAX_INTPOINTS], fz[FEElement::MAX_INTPOINTS];
			for (size_t n = 0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = se.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mpoint->ExtractData<FEElasticMaterialPoint>();
				mat3d m = pt.m_Q;

				// grab the first column as the fiber orientation
				fx[n] = m[0][0];
				fy[n] = m[1][0];
				fz[n] = m[2][0];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			double gx[FEElement::MAX_NODES], gy[FEElement::MAX_NODES], gz[FEElement::MAX_NODES];
			se.project_to_nodes(fx, gx);
			se.project_to_nodes(fy, gy);
			se.project_to_nodes(fz, gz);
			

			for (size_t k = 0; k < neln; k++)
			{
				vec3d v = vec3d(gx[k], gy[k], gz[k]);
				if (m_cultureParams.m_bzfibflat) v.z *= 0.25;
				// normalize the vector
				v.unit();

				// assign the node
				m_pangio->m_fe_node_data[sd.Node(se.m_node[k]).GetID()].m_collfib0 = v;
				m_pangio->m_fe_node_data[sd.Node(se.m_node[k]).GetID()].m_collfib = v;
			}
			//consider setting the fiber orientation of the element here
		}, matls);

		break;
	default:
		//invalid option
		assert(false);
	}
	return true;
}
void FEAngioMaterial::CreateSprouts(double scale)
{
	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_cult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)
	{
		Segment::TIP& tip = *(*tip_it);
		if (tip.bactive)
		{
			AddSprout(tip);
		}
	}
}

void FEAngioMaterial::UpdateSprouts(double scale)
{
	ClearSprouts();
	

	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_cult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)		// Iterate through each segment in the model...
	{
		const Segment::TIP& tip = *(*tip_it);
		assert(tip.bactive);
		assert(tip.pt.ndomain != nullptr);
		assert(tip.pt.elemindex > -1);
		assert(tip.pt.nelem > 0);

		// TODO: What to do with BC==1? Currently, tips that stop growing after hitting boundary
		//       are no longer active. We should still add a sprout for those
		AddSprout(tip);
	}
}

void FEAngioMaterial::UpdateSproutStressScaling()
{
	//TODO: make these user parameters and get better names for these
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	scale = y0 + a / (1 + exp(-(m_pangio->m_time.t - x0) / b));

	return;
}

void FEAngioMaterial::SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp)
{
	FEElasticMaterial::SetLocalCoordinateSystem(el, n, mp);
	FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);

	FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
	FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();

	vessel_elastic.m_Q = pt.m_Q;
	matrix_elastic.m_Q = pt.m_Q;

	FEElasticMaterial* vess_elastic = vessel_material->GetElasticMaterial();
	FEElasticMaterial* mat_elastic = matrix_material->GetElasticMaterial();

	vess_elastic->SetLocalCoordinateSystem(el, n, *angioPt->vessPt);
	mat_elastic->SetLocalCoordinateSystem(el, n, *angioPt->matPt);
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.name(), "sprout") == 0)
	{
		m_suser.push_back(m_s);
	}
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::ClearSprouts()
{
	m_spr.clear();
	sprouts.clear();
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t, FEDomain * domain, int elemindex)
{
	assert(domain != nullptr);
	assert(elemindex != -1);
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(domain);

	vec3d dir = r;
	double pos[3];
	pos[0] = r.x;
	pos[1] = r.y;
	pos[2] = r.z;

	SPROUT s(dir, &dom->Element(elemindex), pos, this);
	assert(s.pel);

	m_spr.push_back(s);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}
//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t, FEDomain * domain)
{
	assert(domain != nullptr);
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESolidDomain * dom = dynamic_cast<FESolidDomain *>(domain);
	double local[3];
	FEElement * el = dom->FindElement(r,local);
	vec3d dir = r;
	double pos[3];
	pos[0] = r.x;
	pos[1] = r.y;
	pos[2] = r.z;

	SPROUT s(dir, el, pos, this);
	assert(s.pel);

	m_spr.push_back(s);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}

void FEAngioMaterial::AddSprout(const Segment::TIP & tip)
{
	double pos[3];
	pos[0] = tip.pt.q.x;
	pos[1] = tip.pt.q.y;
	pos[2] = tip.pt.q.z;

	SPROUT s(tip.u, &tip.pt.ndomain->ElementRef(tip.pt.elemindex), pos, this);

	assert(s.pel);
	m_spr.push_back(s);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}

//-----------------------------------------------------------------------------
vec3d FEAngioMaterial::CurrentPosition(FEElement* pe, double r, double s, double t) const
{
	double arr[FEElement::MAX_NODES];
	FESolidElement * se = dynamic_cast<FESolidElement*>(pe);
	FEMesh * mesh = m_pangio->GetMesh();
	vec3d rc(0,0,0);

	assert(se);
	se->shape_fnc(arr, r, s, t);
	for (int j = 0; j < se->Nodes(); j++)
	{
		rc += mesh->Node(se->m_node[j]).m_rt* arr[j];
	}
	return rc;
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
	FEDomain * d = angioPt.m_pt.ndomain;

	vec3d y;
	assert(angioPt.m_pt.elemindex >= 0);
	y = CurrentPosition(&d->ElementRef(angioPt.m_pt.elemindex), angioPt.m_pt.q.x, angioPt.m_pt.q.x, angioPt.m_pt.q.x);
		
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

			//double p = sprout_s_mag*exp(-sprout_s_range*l);

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
			SPROUT stemp(vec3d(), &angioPt.m_pt.ndomain->Element(angioPt.m_pt.elemindex), local, this);
			temp.push_back(stemp);
			std::pair<size_t, std::vector<SPROUT> *> dim = std::pair<size_t, std::vector<SPROUT> * >(0, &temp);
			std::vector<std::pair<size_t, std::vector<SPROUT> *>> nst = sprouts.within(dim, m_cultureParams.stress_radius * m_cultureParams.stress_radius);
			for (size_t i = 0; i<nst.size(); ++i)
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

void FEAngioMaterial::Update()
{
	m_cult->Update();
}
void FEAngioMaterial::UpdateECM()
{
	ecm_initializer->updateECMdensity(this);
}
bool FEAngioMaterial::Overwrite() const
{
	return ecm_initializer->overwrite();
}
bool FEAngioMaterial::InitCulture()
{
	return m_cult->Init();
}
//this function accumulates the the anistropy and ecm_density, n_tag is incremented to be used to take the average
bool FEAngioMaterial::InitECMDensity(FEAngio * angio)
{
	ecm_initializer->seedECMDensity(this);
	return true;
}
void ECMInitializerConstant::seedECMDensity(FEAngioMaterial * mat)
{
	std::vector<int> matls;
	matls.push_back(mat->GetID());
	mat->m_pangio->ForEachNode([this, mat](FENode & node)
	{
		mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 = mat->m_cultureParams.m_matrix_density;
		mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den  = mat->m_cultureParams.m_matrix_density;
		mat->m_pangio->m_fe_node_data[node.GetID()].m_da = mat->GetAnisotropy();
	}, matls);
}
void ECMInitializer::updateECMdensity(FEAngioMaterial * mat)
{
	std::vector<int> matls;
	matls.push_back(mat->GetID());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	mat->m_pangio->ForEachElement([mesh, mat](FESolidElement & elem, FESolidDomain & d)
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
			nnum = mesh->Node(nnum).GetID();
			// get the ecm density and collagen fiber
			double ecm_den = mat->m_pangio->m_fe_node_data[nnum].m_ecm_den0;
			vec3d coll_fib = mat->m_pangio->m_fe_node_data[nnum].m_collfib0;

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
			mat->m_pangio->m_fe_node_data[nnum].m_collfib += coll_fib;
			mat->m_pangio->m_fe_node_data[nnum].m_ecm_den += ecm_den;


			// increment counter
			mat->m_pangio->m_fe_node_data[nnum].m_ntag++;
		}
	}, matls);
}

void ECMInitializerSpecified::seedECMDensity(FEAngioMaterial * mat)
{
	std::vector<int> matls;
	matls.push_back(mat->GetID());
	FEMesh * mesh = mat->m_pangio->GetMesh();
	mat->m_pangio->ForEachElement([this, mat, mesh](FESolidElement & se, FESolidDomain & sd)
	{
		double den[FEElement::MAX_INTPOINTS];
		double anis[FEElement::MAX_INTPOINTS];
		double pden[FEElement::MAX_INTPOINTS];
		double panis[FEElement::MAX_INTPOINTS];
		for (int n = 0; n<se.GaussPoints(); ++n)
		{
			// generate a coordinate transformation at this integration point
			FEMaterialPoint* mpoint = se.GetMaterialPoint(n);
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
			den[n] = angioPt->m_D;
			anis[n] = angioPt->m_DA;
		}
		se.project_to_nodes(den, pden);
		se.project_to_nodes(anis, panis);
		for (int k = 0; k < se.Nodes(); k++)
		{
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ecm_den0 += pden[k];
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ntag++;
			mat->m_pangio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_da += panis[k];
		}
	}, matls);
}

void ECMInitializerNoOverwrite::seedECMDensity(FEAngioMaterial * mat)
{
	std::vector<int> matls;
	matls.push_back(mat->GetID());
	mat->m_pangio->ForEachNode([this, mat](FENode & node)
	{
		if (mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 == 0.0)
		{
			mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den0 = mat->m_cultureParams.m_matrix_density;
			mat->m_pangio->m_fe_node_data[node.GetID()].m_ecm_den  = mat->m_cultureParams.m_matrix_density;
			mat->m_pangio->m_fe_node_data[node.GetID()].m_da = mat->GetAnisotropy();
		}
	}, matls);
}



double FEAngioMaterial::GetAnisotropy() const
{
	return m_cultureParams.GetWeightInterpolation(1.0);
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
		mat3ds activeStress = AngioStress(*angioPt);
		vessel_elastic.m_s = vessel_material->Stress(*angioPt->vessPt);
		matrix_elastic.m_s = matrix_material->Stress(*angioPt->matPt);
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
		s = angioPt->vessel_weight*vessel_material->Tangent(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->Tangent(*angioPt->matPt);
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
		sed = angioPt->vessel_weight*vessel_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->vessPt) + angioPt->matrix_weight*matrix_material->GetElasticMaterial()->StrainEnergyDensity(*angioPt->matPt);
	}
	return sed;
}
//=============================================================================
BEGIN_PARAMETER_LIST(FEPressureMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_p, FE_PARAM_DOUBLE, "p");
END_PARAMETER_LIST();

mat3ds FEPressureMaterial::Stress(FEMaterialPoint& pt)
{
	mat3dd I(1.0);
	return I*m_p;
}

tens4ds FEPressureMaterial::Tangent(FEMaterialPoint& pt)
{
	return tens4ds(0.0);
}

//============================================================================

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


///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - MirrorSym
//      Calculate force due to mirrored vessels at a particular material point at position x
///////////////////////////////////////////////////////////////////////

void FEAngioMaterial::MirrorSym(vec3d y, mat3ds &s, SPROUT sp, double den_scale)
{
	sym.x = m_cultureParams.m_symmetry_plane.x; sym.y = m_cultureParams.m_symmetry_plane.y; sym.z = m_cultureParams.m_symmetry_plane.z;										// Set the position of the symmetry planes
	mat3ds ssym; ssym.zero();
	vec3d sprout_vect;														// Sprout direction vector
	vec3d sym_v;															// Symmetry vector

	for (int i = 0; i < 7; i++){											// For each of the possible symmetry planes
		if (sym_planes[i] == 1){												// If that symmetry plane is turned on...
			sym_v.x = sym_vects[i][0]; sym_v.y = sym_vects[i][1]; sym_v.z = sym_vects[i][2];	// Obtain the symmetry vector

			vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);
			vec3d r = y - x;																	// Draw the vector r from the sprout location to the position of the material point
			
			//r.x += 2*sym_v.x*(sym.x - x.x); r.y += 2*sym_v.y*(sym.y - x.y); r.z += 2*sym_v.z*(sym.z - x.z);		// Find r for the mirrored vessel sprout
			r.x = r.x + sym_v.x*sym.x; r.y = r.y + sym_v.y*sym.y; r.z = r.z + sym_v.z*sym.z;
			double l = r.unit();													// Find the length of r
			
			sprout_vect.x = sp.sprout.x; sprout_vect.y = sp.sprout.y; sprout_vect.z = sp.sprout.z;	// Set the sprout direction vector 
			sprout_vect.unit();														// Normalize the sprout direction vector
			
			if (m_cultureParams.sprout_s_width != 0){														// If a directional sprout force is being used...
				switch (i) {
				case 0:
					sprout_vect.x = -sprout_vect.x; break;								// Mirror across x
				case 1:
					sprout_vect.y = -sprout_vect.y; break;								// Mirror across y
				case 2:
					sprout_vect.z = -sprout_vect.z; break;								// Mirror across z
				case 3:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; break;	// Mirror across x and y
				case 4:
					sprout_vect.x = -sprout_vect.x; sprout_vect.z = -sprout_vect.z; break;	// Mirror across x and z
				case 5:
					sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;	// Mirror across y and z
				case 6:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;}	// Mirror across x y and z
			}
	
			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double p = den_scale*scale*m_cultureParams.sprout_s_mag*(pow(cos(theta / 2), m_cultureParams.sprout_s_width))*exp(-m_cultureParams.sprout_s_range*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation
			
			if ((p != p) || (r.x != r.x) || (r.y != r.y) || (r.z != r.z)){			// If the mirrored force vector isn't real...
				p = 0.; r.x = 0.; r.y = 0.; r.z = 0.;}									// Set it to zero

			ssym += dyad(r)*p;
		}
	}


	s += ssym;															// Add the symmetry results to the force vector

	return;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEAngioMaterial::CreateMaterialPointData()
{
	return new FEAngioMaterialPoint(new FEElasticMaterialPoint, vessel_material->CreateMaterialPointData(), matrix_material->CreateMaterialPointData()); 
}
