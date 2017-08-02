#include "StdAfx.h"
#include "FEAngioMaterialBase.h"
#include "KDTree/kdtree.h"

#include <FECore/FEElement.h>
#include "FEAngioMaterial.h"

const double PI = 3.141592653589793;
std::vector<double> units3d(3, 1.0);

//-----------------------------------------------------------------------------
vec3d FEAngioMaterialBase::CurrentPosition(FESolidElement * pe, double r, double s, double t) const
{
	double arr[FEElement::MAX_NODES];
	FEMesh * mesh = m_pangio->GetMesh();
	vec3d rc(0, 0, 0);

	assert(pe);
	pe->shape_fnc(arr, r, s, t);
	for (int j = 0; j < pe->Nodes(); j++)
	{
		rc += mesh->Node(pe->m_node[j]).m_rt* arr[j];
	}
	return rc;
}

std::vector<double> access_sprout(std::pair<size_t, std::vector<FEAngioMaterialBase::SPROUT> *> p)
{
	std::vector<double> rv;
	FEAngioMaterialBase::SPROUT * spr = &((*p.second)[p.first]);
	vec3d cpos = spr->mat0->CurrentPosition(spr->pel, spr->r[0], spr->r[1], spr->r[2]);
	rv.emplace_back(cpos.x);
	rv.emplace_back(cpos.y);
	rv.emplace_back(cpos.z);
	return rv;
}

//begin implementation of FEAngioBase
FEAngioMaterialBase::FEAngioMaterialBase() : sprouts(access_sprout, ndim_distance, ndim_distance_to_plane, units3d)
{}

bool FEAngioMaterialBase::FindGridPoint(const vec3d & r, GridPoint & p) const
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

bool FEAngioMaterialBase::FindGridPoint(const vec3d & r, FESolidDomain * domain, int elemindex, GridPoint & p) const
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

FEAngioMaterialBase::SPROUT::SPROUT(const vec3d & dir, FESolidElement * el, double * local, FEAngioMaterialBase * m0, FEElasticMaterial * m1) : sprout(dir), pel(el), mat0(m0), mat1(m1)
{
	r[0] = local[0];
	r[1] = local[1];
	r[2] = local[2];
}

void FEAngioMaterialBase::CreateSprouts(double scale, FEElasticMaterial* emat)
{
	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_cult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)
	{
		Segment::TIP& tip = *(*tip_it);
		if (tip.bactive)
		{
			AddSprout(tip,emat);
		}
	}
}


void FEAngioMaterialBase::UpdateFiberManager()
{
	fiber_manager->Update();
}

void FEAngioMaterialBase::UpdateSprouts(double scale, FEElasticMaterial* emat)
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
		AddSprout(tip,emat);
	}
}

void FEAngioMaterialBase::AddSprout(const Segment::TIP & tip, FEElasticMaterial* emat)
{
	double pos[3];
	pos[0] = tip.pt.q.x;
	pos[1] = tip.pt.q.y;
	pos[2] = tip.pt.q.z;

	m_spr.emplace_back(tip.u, &tip.pt.ndomain->Element(tip.pt.elemindex), pos, this, emat);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}

//-----------------------------------------------------------------------------
void FEAngioMaterialBase::AddSprout(const vec3d& r, const vec3d& t, FEDomain * domain, FEElasticMaterial* emat)
{
	assert(domain != nullptr);
	FESolidDomain * dom = dynamic_cast<FESolidDomain *>(domain);
	double local[3];
	FESolidElement * el = dom->FindElement(r, local);
	vec3d dir = r;
	double pos[3];
	pos[0] = r.x;
	pos[1] = r.y;
	pos[2] = r.z;

	m_spr.emplace_back(dir, el, pos, this, emat);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}

//-----------------------------------------------------------------------------
void FEAngioMaterialBase::ClearSprouts()
{
	m_spr.clear();
	sprouts.clear();
}

//-----------------------------------------------------------------------------
void FEAngioMaterialBase::AddSprout(const vec3d& r, const vec3d& t, FESolidDomain * domain, int elemindex, FEElasticMaterial* emat)
{
	assert(domain != nullptr);
	assert(elemindex != -1);

	vec3d dir = r;
	double pos[3];
	pos[0] = r.x;
	pos[1] = r.y;
	pos[2] = r.z;

	m_spr.emplace_back(dir, &domain->Element(elemindex), pos, this, emat);
	sprouts.insert(std::pair<size_t, std::vector<SPROUT> *>(m_spr.size() - 1, &m_spr));
}

void FEAngioMaterialBase::AdjustMeshStiffness(FEMaterial* mat)
{
	FEMesh & mesh = m_pangio->m_fem->GetMesh();
	if (m_cultureParams.m_composite_material == 0)													// If a composite consitutive model isn't being used, exit
		return;

	int elem_num = 0;													// Element number
	vec3d vess_vect;													// Vessel vector
	std::vector<int> matls;
	matls.emplace_back(mat->GetID());
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1 / static_cast<double>(Nsub);									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

																		//Zero the element items needed
																		//break even on core in field model
	m_pangio->ForEachElementPar([&](FESolidElement & se, FESolidDomain & d)
	{
		int elemnum = se.GetID();
		m_pangio->m_fe_element_data[elemnum].alpha = 0.0;
	}, matls);

	const SegmentList& seg_list = m_cult->GetSegmentList();
	for (ConstSegIter frag_it = seg_list.begin(); frag_it != seg_list.end(); ++frag_it)		// For each segment...
	{
		Segment subunit;												// Segment subdivision placeholder

		const Segment& seg = (*frag_it);												// Obtain the segment

		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			assert(seg.length() > 0.);
			if (k == 1) {													// If this is the first subdivision, find the origin of the segment
				subunit.tip(0).pt = seg.tip_c(0).pt;
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
				m_pangio->IsInsideHex8(&static_cast<FESolidElement&>(subunit.tip(1).pt.ndomain->ElementRef(subunit.tip(1).pt.elemindex)), mid, &mesh, intp))
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
			}

			// Set the origin of the next subdivision to the end of the current one
			subunit.tip(0).pt.r = subunit.tip(1).pos();
		}
	}

	// Volume fraction for the composite material model



	//verified good on core in field model
	m_pangio->ForEachElementPar([&](FESolidElement & e, FESolidDomain & d)
	{
		assert(std::find(domainptrs.begin(), domainptrs.end(), &d) != domainptrs.end());
		vec3d e1; vec3d e2; vec3d e3;						// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
		int nint = e.GaussPoints();										// Obtain the number of gauss points
		int num_elem = e.GetID();

		FEAngioElementData& eg = m_pangio->m_fe_element_data[num_elem]; // Obtain the element from the domain	
		double alpha = eg.alpha;										// Obtain alpha from the grid element
		for (int n = 0; n< e.Nodes(); n++)
		{
			int id = e.m_node[n];
			id = mesh.Node(id).GetID();
			m_pangio->m_fe_node_data[id].alpha = alpha;
		}

		for (int n = 0; n < nint; ++n)										// For each gauss point...
		{
			FEMaterialPoint& mp = *(e.GetMaterialPoint(n));
			FEAngioMaterialPoint* pt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp); // get the mixture material point
			pt->vessel_weight = alpha;
			pt->matrix_weight = 1.0 - alpha;
		}
	}, matls);

	return;
}

///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - MirrorSym
//      Calculate force due to mirrored vessels at a particular material point at position x
///////////////////////////////////////////////////////////////////////

void FEAngioMaterialBase::MirrorSym(vec3d y, mat3ds &s, SPROUT sp, double den_scale)
{
	sym.x = m_cultureParams.m_symmetry_plane.x; sym.y = m_cultureParams.m_symmetry_plane.y; sym.z = m_cultureParams.m_symmetry_plane.z;										// Set the position of the symmetry planes
	mat3ds ssym; ssym.zero();
	vec3d sprout_vect;														// Sprout direction vector
	vec3d sym_v;															// Symmetry vector

	for (int i = 0; i < 7; i++) {											// For each of the possible symmetry planes
		if (sym_planes[i] == 1) {												// If that symmetry plane is turned on...
			sym_v.x = sym_vects[i][0]; sym_v.y = sym_vects[i][1]; sym_v.z = sym_vects[i][2];	// Obtain the symmetry vector

			vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);
			vec3d r = y - x;																	// Draw the vector r from the sprout location to the position of the material point

																								//r.x += 2*sym_v.x*(sym.x - x.x); r.y += 2*sym_v.y*(sym.y - x.y); r.z += 2*sym_v.z*(sym.z - x.z);		// Find r for the mirrored vessel sprout
			r.x = r.x + sym_v.x*sym.x; r.y = r.y + sym_v.y*sym.y; r.z = r.z + sym_v.z*sym.z;
			double l = r.unit();													// Find the length of r

			sprout_vect.x = sp.sprout.x; sprout_vect.y = sp.sprout.y; sprout_vect.z = sp.sprout.z;	// Set the sprout direction vector 
			sprout_vect.unit();														// Normalize the sprout direction vector

			if (m_cultureParams.sprout_s_width != 0) {														// If a directional sprout force is being used...
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
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;
				}	// Mirror across x y and z
			}

			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double p = den_scale*scale*m_cultureParams.sprout_s_mag*(pow(cos(theta / 2), m_cultureParams.sprout_s_width))*exp(-m_cultureParams.sprout_s_range*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

			if ((p != p) || (r.x != r.x) || (r.y != r.y) || (r.z != r.z)) {			// If the mirrored force vector isn't real...
				p = 0.; r.x = 0.; r.y = 0.; r.z = 0.;
			}									// Set it to zero

			ssym += dyad(r)*p;
		}
	}


	s += ssym;															// Add the symmetry results to the force vector

	return;
}
