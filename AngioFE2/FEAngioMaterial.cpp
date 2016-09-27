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
	ADD_PARAMETER2(m_cultureParams.m_anastomosis_distance, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "anastomosis_distance");
	ADD_PARAMETER2(m_cultureParams.m_branch_chance, FE_PARAM_DOUBLE, FE_RANGE_RIGHT_OPEN(0,1), "branch_chance");
	ADD_PARAMETER(m_cultureParams.m_branch, FE_PARAM_BOOL, "branch");
	ADD_PARAMETER(m_cultureParams.m_anastomosis, FE_PARAM_BOOL, "anastomosis");
	ADD_PARAMETER2(m_cultureParams.m_vessel_width, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "vessel_width");
	ADD_PARAMETER(m_cultureParams.m_boundary_condition_type, FE_PARAM_STRING, "boundary_condition_type");

	ADD_PARAMETER2(m_cultureParams.m_matrix_condition, FE_PARAM_INT, FE_RANGE_CLOSED(0,4), "matrix_condition");
	ADD_PARAMETER2(m_cultureParams.m_matrix_density, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "matrix_density");

	ADD_PARAMETER(m_cultureParams.m_symmetry_plane, FE_PARAM_VEC3D, "symmetryplane");
	//uncategorized variables are incomplete
	ADD_PARAMETER(m_cultureParams.m_composite_material, FE_PARAM_INT, "composite_material");
	ADD_PARAMETER2(m_cultureParams.m_sprout_force, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "sprout_force");
	ADD_PARAMETER2(m_cultureParams.m_number_fragments, FE_PARAM_INT, FE_RANGE_GREATER_OR_EQUAL(0), "number_fragments" );
	ADD_PARAMETER(m_cultureParams.vessel_orient_weights, FE_PARAM_VEC3D, "weights");
	ADD_PARAMETER(m_cultureParams.m_seed, FE_PARAM_INT, "seed");
	
END_PARAMETER_LIST();
//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	scale = 1.0;
	
	m_pangio = nullptr;

	m_cult = nullptr;

	AddProperty(&vessel_material, "vessel");
	AddProperty(&matrix_material , "matrix");
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


	if(matrix_material->Init() == false) return false;

	if(vessel_material->Init() == false) return false;

	if (FEElasticMaterial::Init() == false) return false;


	//culture must be initialized here  so pangio is defined
	assert(m_pangio);
	m_cult = new Culture(*m_pangio, this, &m_cultureParams);

	// add the user sprouts
	std::vector<int> matls;
	matls.push_back(this->GetID());
	FEMesh& mesh = GetFEModel()->GetMesh();
	mesh.DomainListFromMaterial(matls, domains);
	for (auto i = 0; i < domains.size(); i++)
	{
		domainptrs.push_back(&mesh.Domain(domains[i]));
	}
	for (unsigned int i=0; i<m_suser.size(); ++i)
	{
		if (domains.size())
			AddSprout(m_suser[i], vec3d(0,0,0), &mesh.Domain(domains[0]));
		//TODO: sprouts probably need distributed among the domains of the material
	}
	m_suser.clear();

	// initialize material point data
	vec3d x[FEElement::MAX_NODES];
	
	for (int n=0; n<mesh.Domains(); ++n)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(n));
		FEMaterial* pm = dom.GetMaterial();
		FEAngioMaterial* pam;
		if(strcmp(pm->GetTypeStr(), "angio")==0)
		{
			pam = dynamic_cast<FEAngioMaterial*>(pm);
		}
		else
		{
			pam = dynamic_cast<FEAngioMaterial*>(pm->FindComponentByType("angio"));
		}
		if (pam == this)
		{
			// loop over all elements
			int NE = dom.Elements();
			for (int i=0; i<NE; ++i)
			{
				// get the next element
				FEElement& el = dom.Element(i);
				int neln = el.Nodes();
				
				// get the nodal coordinates
				//TODO: understand why m_node is used, these are the adjacent nodes
				for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

				// loop over all integration points
				int nint = el.GaussPoints();
				for (int j=0; j<nint; ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);
					FEAngioMaterialPoint* pt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
					if(pt)
					{
						vec3d r = el.Evaluate(x, j);

						// calculate the GridPoint data for this point.
						//TODO: check elastic material for integration point coordinates
						if (m_pangio->FindGridPoint(r, &dom,i, pt->m_pt) == false)
						{
							return false;
						}
					}
				}
			}
		}
	}

	return true;
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
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length() > 0.){											// If it's a +1 segment...
					subunit.tip(0).pt = seg.tip(0).pt;
				}
				else{															// If it's a -1 segment...
					subunit.tip(1).pt.r = seg.tip(1).pos();
					//					subunit.m_length = -1.;
					assert(false);
				}
			}

			// Calculate the subdivision
			vec3d v = seg.uvect();
			if (seg.length() > 0.){
				subunit.tip(1).pt.r = subunit.tip(0).pos() + v*(sub_scale*seg.length());
			}
			else{														// If it's a -1 segment...
				subunit.tip(0).pt.r = subunit.tip(1).pos() + v*(sub_scale*seg.length());
				assert(false);
			}

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
			if (seg.length() > 0.){										// If it's a +1 segment...
				vess_vect = subunit.tip(1).pos() - subunit.tip(0).pos();
			}
			else{														// If it's a -1 segment...
				vess_vect = subunit.tip(0).pos() - subunit.tip(1).pos();
			}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
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
			if (seg.length() > 0.){
				subunit.tip(0).pt.r = subunit.tip(1).pos();
			}
			else{
				subunit.tip(1).pt.r = subunit.tip(0).pos();
			}
		}
	}

	// Volume fraction for the composite material model




	m_pangio->ForEachElement([this, &mesh](FESolidElement & e, FESolidDomain & d)
	{
		vec3d e1; vec3d e2; vec3d e3;						// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
		double alpha = 0.;									// Obtain the element from the domain
		int nint = e.GaussPoints();										// Obtain the number of gauss points
		int num_elem = e.GetID();

		FEAngioElementData& eg = m_pangio->m_fe_element_data[num_elem];
		alpha = eg.alpha;											// Obtain alpha from the grid element
		for (int n = 0; n< e.Nodes(); n++)
		{
			//hack
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

bool FEAngioMaterial::InitCollagenFibers()
{
	std::vector<int> matls;
	matls.push_back(GetID());
	//remove later
	FEMesh * mesh = m_pangio->GetMesh();
	std::vector<vec3d> colfibs;
	colfibs.reserve(mesh->Nodes());

	switch (m_cultureParams.m_matrix_condition)
	{
	case 0: // random orientation
		//TODO: loop replaced for consistent answers
		//TODO: refactor once code is working
		
		for (int i = 0; i < mesh->Nodes(); i++)
		{
			vec3d v = vrand();
			if (m_cultureParams.m_bzfibflat) v.z *= 0.25;

			// normalize the vector
			v.unit();
			colfibs.push_back(v);
		}

		for (int i = 0; i < mesh->Nodes();i++)

		//m_pangio->ForEachNode([&](FENode & node)
		{
			

			// assign the node
			FENode & node = mesh->Node(i);
			
			//m_pangio->m_fe_node_data[node.GetID()].m_collfib0 = v;
			//m_pangio->m_fe_node_data[node.GetID()].m_collfib = v;
			//hack
			vec3d v = colfibs[i];
			m_pangio->m_fe_node_data[node.GetID()].m_collfib0 = v;
			m_pangio->m_fe_node_data[node.GetID()].m_collfib = v;
		}//, matls);
		
		/*
		m_pangio->ForEachNode([&](FENode & node)
		{
			if (m_pangio->m_fe_node_data[node.GetID()].m_collfib0.norm() < 0.5)
			{
				vec3d v = vrand();
				if (m_cultureParams.m_bzfibflat) v.z *= 0.25;

				// normalize the vector
				v.unit();

				m_pangio->m_fe_node_data[node.GetID()].m_collfib0 = v;
				m_pangio->m_fe_node_data[node.GetID()].m_collfib = v;
			}
			
		}, matls);
		*/
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
	double magnitude = scale*m_cultureParams.m_sprout_force;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_pangio->m_time.t == 0.0)
		magnitude = (1.0 / 4.0)*0.001*scale;
	else if (m_pangio->m_time.t < 4.0)
		magnitude = (1.0 / 4.0)*m_pangio->m_time.t*scale;

	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_cult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)
	{
		Segment::TIP& tip = *(*tip_it);
		if (tip.bactive)
		{
			// get the tip
			const vec3d& r = tip.pos();

			// get the directional unit vector of the tip
			const vec3d& u = tip.u;

			AddSprout(r, u, tip.pt.ndomain);
			tip.bdyf_id = Sprouts() - 1;
		}
	}
}

void FEAngioMaterial::UpdateSprouts(double scale)
{
	double magnitude = scale*m_cultureParams.m_sprout_force;								// Magnitude of the sprout force
	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_pangio->m_time.t == 0.0)
		magnitude = (1.0 / 4.0)*0.001*scale;
	else if (m_pangio->m_time.t < 4.0)
		magnitude = (1.0 / 4.0)*m_pangio->m_time.t*scale;

	ClearSprouts();

	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_cult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)		// Iterate through each segment in the model...
	{
		const Segment::TIP& tip = *(*tip_it);
		assert(tip.bactive);
		assert(tip.bdyf_id >= 0);

		// TODO: What to do with BC==1? Currently, tips that stop growing after hitting boundary
		//       are no longer active. We should still add a sprout for those
		AddSprout(tip.pos(), tip.u, tip.pt.ndomain);
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
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t, FEDomain * domain, int elemindex)
{
	assert(domain != nullptr);
	assert(elemindex != -1);
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESolidDomain* dom = dynamic_cast<FESolidDomain*>(domain);

	SPROUT s;
	s.pel = &(dom->ElementRef(elemindex));
	s.sprout = t;

	assert(s.pel);

	m_spr.push_back(s);
}
//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t, FEDomain * domain)
{
	assert(domain != nullptr);
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESolidDomain * dom = dynamic_cast<FESolidDomain *>(domain);

	SPROUT s;
	s.pel = dom->FindElement(r,s.r);
	s.sprout = t;

	assert(s.pel);

	m_spr.push_back(s);
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
		

	//#pragma omp parallel for shared(s)
	for (int i=0; i<NS; ++i)
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

		if (sym_on == true)															// If symmetry is turned on, apply symmetry
			MirrorSym(y, si, sp, den_scale);

//#pragma omp critical
		s += si;
	}
	return s;
}

void FEAngioMaterial::Grow(SimulationTime& time)
{
	m_cult->Grow(time);
}
void FEAngioMaterial::Update()
{
	m_cult->Update();
}
bool FEAngioMaterial::InitCulture()
{
	return m_cult->Init();
}
//this function accumulates the the anistropy and ecm_density, n_tag is incremented to be used to take the average
bool FEAngioMaterial::InitECMDensity(FEAngio * angio)
{
	if (m_cultureParams.m_matrix_density == 0.0)
	{
		std::vector<int> matls;
		matls.push_back(GetID());
		FEMesh * mesh = angio->GetMesh();
		angio->ForEachElement([this, angio, mesh](FESolidElement & se, FESolidDomain & sd)
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
				angio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ecm_den0 += pden[k];
				angio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_ntag++;
				angio->m_fe_node_data[mesh->Node(se.m_node[k]).GetID()].m_da += panis[k];
			}
		}, matls);
	}
	else
	{
		std::vector<int> matls;
		matls.push_back(GetID());
		angio->ForEachNode([this, angio](FENode & node)
		{
			angio->m_fe_node_data[node.GetID()].m_ecm_den0 += this->m_cultureParams.m_matrix_density;
			angio->m_fe_node_data[node.GetID()].m_ntag++;
			angio->m_fe_node_data[node.GetID()].m_da = this->GetAnisotropy();
		}, matls);
	}

	return true;
}

double FEAngioMaterial::GetAnisotropy() const
{
	double ans = m_cultureParams.vessel_orient_weights.y / (m_cultureParams.vessel_orient_weights.x + m_cultureParams.vessel_orient_weights.y);
	return ans;
}
void FEAngioMaterial::SetBoundaryCondition() const
{
	switch (m_cultureParams.m_boundary_condition_type[0])
	{
	case 's':
		m_cult->ChangeBC(*m_pangio, BC::STOP);
		break;
	case 'b':
		m_cult->ChangeBC(*m_pangio, BC::BOUNCY);
		break;
	default:
		assert(false);
	}
}

mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
	FEElasticMaterialPoint& elastic_pt = *mp.ExtractData<FEElasticMaterialPoint>();

	mat3ds s;
	s.zero();
	//should always be true but we should check
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
