#include "StdAfx.h"
#include "Culture.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"
#include <random>
#define _USE_MATH_DEFINES

void DirectionalWeights(double da, double dw[2]);

CultureParameters::CultureParameters()
{
	m_boundary_condition_type[0] = 's';
	m_boundary_condition_type[1] = '\0';
}

//-----------------------------------------------------------------------------
Culture::Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp) : m_angio(angio)
{
	m_cultParams = cp;

	m_num_zdead = 0;

	// TODO: I think it would make sense to tie the anastomosis distance 
	//       to the vessel diameter.
    m_total_length = 0.;
	m_num_vessel = 0;
	m_num_branches = 0;		// Initialize branching counter

	// Initialize counters
	m_num_anastom = 0;
	bc = nullptr;
	ChangeBC(angio, BC::STOP);
	assert(bc);
	m_nsegs = 0;			// Initialize segment counter
	m_pmat = matl;
}

void Culture::ChangeBC(FEAngio & angio, int bcset)
{
	if (bc)
	{
		delete bc;
		bc = nullptr;
	}
		
	switch (bcset)
	{
	case BC::STOP:
		bc = new StopBC(angio, this);
		break;
	case BC::BOUNCY:
		bc = new BouncyBC(angio, this);
		break;
	}
}

//-----------------------------------------------------------------------------
Culture::~Culture()
{
	delete bc;
	if (fseeder)
		delete fseeder;
}

//-----------------------------------------------------------------------------
// Initialize the culture
bool Culture::Init()
{
	//check culture parameters for valid values here


	double d = m_cultParams->m_y0 + m_cultParams->m_culture_a / (1.0 + exp(m_cultParams->m_x0 / m_cultParams->m_culture_b)); // Initial value of growth curve (t = 0)
	m_vess_length = d;

	// If branching is turned off, we set the branching chance to zero
	if (m_cultParams->m_branch == false) m_cultParams->m_branch_chance = 0.0;

	// make sure the initial length is initialized 
	if (m_cultParams->m_initial_vessel_length <= 0.0) m_cultParams->m_initial_vessel_length = m_vess_length;

	//intialize the Fragment Seeder
	switch (m_cultParams->fragment_seeder)
	{
	case 0:
		fseeder = new ClassicFragmentSeeder(m_cultParams, m_angio);
		break;
	case 1:
		fseeder = new MultiDomainFragmentSeeder(m_cultParams, m_angio);
		break;
	default:
		assert(false);
	}

	// do the initial seeding
	return fseeder->SeedFragments(m_angio.CurrentSimTime(), this);
}
//currently all weights are set to 1 but in the future override this and change it to use the segment volume or something else 
void Culture::SetWeights(vector<SegGenItem> & weights, std::vector<FEDomain*> & domains)
{
	FEMesh * mesh = m_angio.GetMesh();
	for (size_t i = 0; i < domains.size(); i++)
	{
		for (auto j = 0; j < domains[i]->Elements(); j++)
		{
			SegGenItem cur;
			cur.weight = 1.0;
			cur.domain = domains[i];
			cur.ielement = j;
			weights.push_back(cur);
		}
	}
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool ClassicFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	for (int i = 0; i < culture_params->m_number_fragments; ++i)
	{
		// Create an initial segment
		Segment seg;
		if (createInitFrag(seg) == false) return false;

		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		// add it to the list
		culture->AddSegment(seg);
	}

	// init vessel counter
	culture->m_num_vessel = culture_params->m_number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	return true;
}

//-----------------------------------------------------------------------------
// Generates an initial fragment that lies inside the grid.
bool ClassicFragmentSeeder::createInitFrag(Segment& seg)
{
	FEMesh * mesh = m_angio.GetMesh();
	// Set seg length to value of growth function at t = 0
	double seg_length = culture_params->m_initial_vessel_length;

	// We create only segments that lie inside the grid.
	// Since the creation of such a segment may fail (i.e. too close to boundary)
	// we loop until we find one.
	const int MAX_TRIES = 10;
	int ntries = 0;
	GridPoint p0, p1;
	do
	{
		// pick a random element
		// The factor 0.999 is because there is a small chance that
		// frand() returns 1 in which case elem_num would be invalid.
		int elem_num = int(0.999*frand()*mesh->Domain(0).Elements());
		//remove these when error checking is complete
		//std::vector<int> validElements = {283,286, 253, 256};
		//int elem_num = validElements[rand() % validElements.size()];
		// generate random natural coordinates
		vec3d q = vrand();

		// set the position of the first tip
		p0 = m_angio.FindGridPoint(&mesh->Domain(0), elem_num, q);

		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = vrand();
		seg_vec.unit();

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		m_angio.m_pmat[0]->FindGridPoint(r1,p1);

		ntries++;
	} while ((p1.nelem == -1) && (ntries < MAX_TRIES));

	if (p1.nelem == -1)  return false;

	// assign the grid points
	seg.tip(0).pt = p0;
	seg.tip(1).pt = p1;

	// update length and unit vector
	seg.Update();

	// make both tips active
	seg.tip(0).bactive = true;
	seg.tip(1).bactive = true;

	// decide if this initial segment is allowed to branch
	if (frand() < culture_params->m_initial_branch_probability) seg.SetFlagOn(Segment::INIT_BRANCH);

	// Mark segment as an initial fragment
	seg.SetFlagOn(Segment::INIT_SPROUT);

	// all good
	return true;
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool MultiDomainFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	std::default_random_engine reg;
	reg.seed(m_angio.m_irseed);
	FEMesh * mesh = m_angio.GetMesh();
	int elementsInDomains = 0;
	for (size_t i = 0; i < culture->m_pmat->domains.size(); i++)
	{
		elementsInDomains += m_angio.GetMesh()->Domain(culture->m_pmat->domains[i]).Elements();
		domains.push_back(&m_angio.GetMesh()->Domain(culture->m_pmat->domains[i]));
	}
		
	std::uniform_int_distribution<int> ddist(0, culture->m_pmat->domains.size() - 1);
	std::uniform_int_distribution<int> edist(0, elementsInDomains -1);
	SegGenItem sgi;

	for (int i = 0; i < culture_params->m_number_fragments; ++i)
	{
		// Create an initial segment
		Segment seg;
		sgi.domain = &mesh->Domain(culture->m_pmat->domains[ddist(reg)]);
		sgi.ielement = edist(reg);
		if (createInitFrag(seg, sgi, culture) == false) return false;

		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		// add it to the list
		culture->AddSegment(seg);
	}

	// init vessel counter
	culture->m_num_vessel = culture_params->m_number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	return true;
}

//-----------------------------------------------------------------------------
// Generates an initial fragment that lies inside the given element
bool MultiDomainFragmentSeeder::createInitFrag(Segment& seg, SegGenItem & item, Culture * culture)
{
	//note tip(1) may not be in the same element as the initial fragment
	// Set seg length to value of growth function at t = 0
	double seg_length = culture_params->m_initial_vessel_length;

	// Since the creation of such a segment may fail (i.e. too close to boundary)
	// we loop until we find one.
	const int MAX_TRIES = 10;
	int ntries = 0;
	GridPoint p0, p1;
	do
	{
		// generate random natural coordinates
		vec3d q = vrand();

		// set the position of the first tip
		p0 = m_angio.FindGridPoint(item.domain, item.ielement, q);

		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = vrand();
		seg_vec.unit();

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		culture->m_pmat->FindGridPoint(r1,p1);
		//need to check the domain is legal
		ntries++;
	} while (((p1.nelem == -1) && (ntries < MAX_TRIES)) || (std::find(domains.begin(), domains.end(), p1.ndomain) == domains.end()));

	if (p1.nelem == -1)  return false;

	// assign the grid points
	seg.tip(0).pt = p0;
	seg.tip(1).pt = p1;

	// update length and unit vector
	seg.Update();

	// make both tips active
	seg.tip(0).bactive = true;
	seg.tip(1).bactive = true;

	// decide if this initial segment is allowed to branch
	if (frand() < culture_params->m_initial_branch_probability) seg.SetFlagOn(Segment::INIT_BRANCH);

	// Mark segment as an initial fragment
	seg.SetFlagOn(Segment::INIT_SPROUT);

	// all good
	return true;
}

//-----------------------------------------------------------------------------
// Vessel elongation is represented by the addition of new line segments at the locations of the active sprouts.
void Culture::Grow(SimulationTime& time)
{
	// Determine length of new segments for this growth step        		
	UpdateNewVesselLength(time);

	// Elongate the active vessel tips
	GrowVessels();
	
	// Branching phase
	if (m_cultParams->m_branch) BranchVessels(time);

	// Anastomosis phase
	if (m_cultParams->m_anastomosis) FuseVessels();

	// Update all active growth tips
	FindActiveTips();
}

//-----------------------------------------------------------------------------
// The new vessel length (i.e. the length of new vessels) is a function
// of time. This function is called at the start of each growth step (Culture::Grow)
// and updates the new vessel length.
void Culture::UpdateNewVesselLength(SimulationTime& time)
{
	double lc = m_cultParams->m_culture_a/(1.+ exp(-(time.t - m_cultParams->m_x0)/m_cultParams->m_culture_b));
	lc -= m_cultParams->m_culture_a/(1. + exp(-(time.t - time.dt-m_cultParams->m_x0)/m_cultParams->m_culture_b));
	
	m_vess_length = lc*m_cultParams->m_length_adjustment;
}

//-----------------------------------------------------------------------------
// Grow phase
void Culture::GrowVessels()
{
	FEMesh* mesh = m_angio.GetMesh();

	for (TipIter it = m_active_tips.begin(); it != m_active_tips.end(); ++it)
	{
		Segment::TIP& tip = *(*it);

		//modify direction based on anisotropy value
		double dw[2];
		double da_value = m_angio.genericProjectToPoint(*tip.pt.ndomain, &tip.pt.ndomain->ElementRef(tip.pt.elemindex), &FEAngioNodeData::m_da, tip.pt.q.x, tip.pt.q.y, tip.pt.q.z);
		DirectionalWeights(da_value, dw);
		//TODO: this overwrite may not be ideal
		m_cultParams->vessel_orient_weights.x = dw[0];
		m_cultParams->vessel_orient_weights.y = dw[1];
		
		// Create new vessel segment at the current tip existing segment
		Segment seg = GrowSegment(tip);

		// Add the segment to the network
		// This will also enforce the boundary conditions
		AddNewSegment(seg);
    }
}

//-----------------------------------------------------------------------------
// Create a new segment at the (active) tip of an existing segment
Segment Culture::GrowSegment(Segment::TIP& tip, bool branch, bool bnew_vessel, vec3d growthDirection)
{

	// Make sure the tip is active
	assert(tip.bactive);

	// calculate the length scale factor based on collagen density
	double den_scale = FindDensityScale(tip.pt);

	// this is the new segment length
	double seg_length = den_scale*m_vess_length;

	// determine the growth direction
	vec3d seg_vec;
	if(growthDirection.x==0 && growthDirection.y==0 && growthDirection.z == 0)
	{
		seg_vec = FindGrowDirection(tip);
	}
	else
	{
		seg_vec = growthDirection;
	}

	// If new segment is a branch we modify the grow direction a bit
	if (branch)
	{
		// TODO: what's the logic here? Why the 0.5 factor?
		//      If the vessel is aligned with the collagen (and the initial fragments are)
		//      then  the new branch will overlap the old segment. 
		vec3d coll_fib = m_angio.CollagenDirection(tip.pt);
		seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
		seg_vec.unit();

		// TODO: Trying something new here: grow perpendicular to segment.
		//vec3d e = vrand(); e.unit();
		//double d = e*seg_vec;
		//assert(d != 0.0);
		//seg_vec = e - seg_vec*d;
		//seg_vec.unit();
	}

	// Create a new segment
	Segment seg;

	// transer seed label
	seg.seed(tip.nseed);

	// assign vessel ID
	if (bnew_vessel)
	{
		seg.vessel(m_num_vessel++);
	}
	else seg.vessel(tip.nvessel);

	// create a new vessel
	// the first point is a copy of the last point of the source segment
	seg.tip(0).pt = tip.pt;

	// position the end point
	seg.tip(1).pt.r = seg.tip(0).pos() + seg_vec*seg_length;

	seg.tip(0).bactive = false;		// Turn off origin tip of new segment
    seg.tip(1).bactive = true;		// Turn on end tip of new segment

	// Turn off previous segment tip
	tip.bactive = false;


	//TODO: check if still needed
	if (branch)
	{
		m_pmat->FindGridPoint(seg.tip(1).pos(),seg.tip(1).pt);
		//should fill in the  fields of the new segment

		//do a check if the new vessel is in an element
	}

	// update length and unit vector
	seg.Update();

	return seg;
}

//-----------------------------------------------------------------------------
// Branching phase of the growth step.
void Culture::BranchVessels(SimulationTime& time)
{
	// Elongate the active vessel tips
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		// decide if this segment can branch
		bool branch = true;
		if (it->GetFlag(Segment::BC_DEAD)) branch = false;	// Make sure segment has not reached a boundary
		if (it->GetFlag(Segment::ANAST  )) branch = false;	// and has not undegone anastimoses.
		if (it->tip(0).bactive || it->tip(1).bactive) branch = false;	// don't branch segments that are actively growing
		int vessel = it->m_nvessel;
		if (branch || it->GetFlag(Segment::INIT_BRANCH))
		{
			// pick an end randomly
			int k = ( frand() < 0.5 ? k = 0 : k = 1);

			// find the density scale factor
			double den_scale = FindDensityScale(it->tip(k).pt);

			// calculate actual branching probability
			// TODO: This formula makes no sense to me.
			double t = time.t;
			double dt = time.dt;
			double bprob = den_scale*dt*m_cultParams->m_branch_chance/t;
			assert(bprob < 1.0);

			// If branching is turned on determine if the segment forms a new branch
			if (frand() < bprob)
			{
				// Turn off the initial branch flag
				it->SetFlagOff(Segment::INIT_BRANCH);
				BranchSegment(it->tip(k));

				// Increase the number of branches.
				m_num_branches++;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Branching is modeled as a random process, determine is the segment passed to this function forms a branch
void Culture::BranchSegment(Segment::TIP& tip)
{
	// we must reactive the tip
	// (it will be deactivated by GrowSegment)
	tip.bactive = true;

	// Create the new vessel segment (with branch flag and new_vessel flag to true)
	Segment seg = GrowSegment(tip, true, true);

	// Add it to the culture
	AddNewSegment(seg);
}

//-----------------------------------------------------------------------------
// Create a new sprout force component for a newly formed branch
// TODO: What if both tips are active?
void Culture::CreateBranchingForce(Segment& seg)
{
	m_angio.total_bdyf = 0;																// Obtain the total number of sprouts
	if (m_pmat) m_angio.total_bdyf = m_pmat->Sprouts();
	else m_angio.total_bdyf = m_angio.m_pbf->Sprouts();

	vec3d tip, sprout_vect;
	vec3d seg_vec = seg.uvect();
	Segment::TIP stip;
	
	if (seg.tip(0).bactive)
	{
		stip = seg.tip(0);
		tip = seg.tip(0).pos();															// Obtain the position of the new tip
		sprout_vect = -seg_vec;	// notice negative sign										
	}
	else if (seg.tip(1).bactive)
	{
		stip = seg.tip(1);
		tip = seg.tip(1).pos();
		sprout_vect = seg_vec;
	}
	else
	{
		//what to do when neither tip is active?
		//this is ussually reached when the branch has grown out of bounds
		stip = seg.tip(1);
		tip = seg.tip(1).pos();
		sprout_vect = seg_vec;
	}

	if (m_pmat)
	{
		m_pmat->AddSprout(tip, sprout_vect, stip.pt.ndomain, stip.pt.elemindex);
		m_angio.total_bdyf = m_pmat->Sprouts();
	}
	else
	{
		m_angio.m_pbf->AddSprout(tip, sprout_vect);						// Add the new sprout component to the sprout force field
		m_angio.total_bdyf = m_angio.m_pbf->Sprouts();												// Update the total number of sprouts
	}
}

//-----------------------------------------------------------------------------
// Anastimoses phase.
// TODO: This only implements fusing at tips. Maybe we should extend this to do line-line intersections.
void Culture::FuseVessels()
{
	// loop over all segments
    for (SegIter it1 = m_frag.begin(); it1 != m_frag.end(); ++it1)
	{
		// Make sure the vessels has not fused yet
		if (it1->GetFlag(Segment::ANAST) == false)
		{
			// loop over tips
			for (int k1=0; k1<2; ++k1)
			{
				// Make sure the tip is active
				if (it1->tip(k1).bactive)
				{
					for (SegIter it2 = m_frag.begin(); it2 != m_frag.end(); ++it2)
					{
						// make sure neither segments sprout from the same initial fragment
						// TODO: why is this not allowed? I think this might be to prevent
						//       anastimoses after branching since the branched vessels can be small. Find a better way!
						if (it1->seed() != it2->seed())
						{
							double dist0 = (it1->tip(k1).pos() - it2->tip(0).pos()).norm();
							double dist1 = (it1->tip(k1).pos() - it2->tip(1).pos()).norm();

							// pick the closest tip
							double dist = dist0;
							int k2 = 0;
							if (dist1 < dist0) { k2 = 1; dist = dist1; }

							// see if we can fuse
							if (dist < m_cultParams->m_anastomosis_distance)
							{
								// create a segment between the two segments to complete the anastomosis
								Segment seg = ConnectSegment(*it1, *it2, k1, k2);

								// mark segments
								seg.SetFlagOn(Segment::ANAST);
								it1->SetFlagOn(Segment::ANAST);
								it2->SetFlagOn(Segment::ANAST);

								if (bc->ChangeOfMaterial(seg))
									break;
								// add it to the list
								if (seg.length() > 0.0)
									AddSegment(seg);

								// deactivate tips
								it1->tip(k1).bactive = false;
								it2->tip(k2).bactive = false;	// for end-to-end anastimoses.
							
								// increment counter
								++m_num_anastom;

								// break it2 loop since tip(k1) is no longer active
								break;
							}
						}
					}
				}
			}
		}
	} 
}

//-----------------------------------------------------------------------------
// Add a new segment to the culture.
// This will apply BCs to the new segment and may result in 
// the addition of several new segments. 
// (it is assumed that the other tip was form at the end of 
// another segment and is valid)
// we always assume that tip(0) is the "old" tip and tip(1)
// is the new tip
void Culture::AddNewSegment(Segment& seg)
{
	//adding zero length segments shoudl be avoided
	assert(seg.length() > 0.0);
	assert(seg.tip(0).pt.nelem >= 0);

	// get the new tip
	Segment::TIP& new_tip = seg.tip(1);
	//assert(new_tip.pt.nelem == -1);
	//init done elsewhere
	assert(new_tip.bactive);

	bc->CheckBC(seg);
}

//-----------------------------------------------------------------------------
// This function adds a segment to the list.
// It also updates the segment counter and assigns Segment::m_nid and sets the time of birth.
// Note that we also push to the front of the list so that we don't corrupt
// any active iterators.
void Culture::AddSegment(Segment& seg)
{
	// let's check a few things
	assert(seg.tip(0).pt.nelem >= 0);
	assert(seg.tip(1).pt.nelem >= 0);

	assert(seg.tip(0).pt.ndomain != nullptr);
	assert(seg.tip(1).pt.ndomain != nullptr);

	assert(seg.tip(0).pt.elemindex >= 0);
	assert(seg.tip(1).pt.elemindex >= 0);

	assert(seg.tip(0).nseed >= 0);
	assert(seg.tip(1).nseed >= 0);

	assert(seg.tip(0).nvessel >= 0);
	assert(seg.tip(1).nvessel >= 0);

	assert(seg.seed() >= 0);
	assert(seg.vessel() >= 0);

	assert(seg.length() > 0.0);

	vec3d rpo = m_angio.Position(seg.tip(0).pt);
	vec3d rp = seg.tip(0).pt.r;
	double dist = distance(rpo.x, rp.x, rpo.y, rp.y, rpo.z, rp.z);
	assert(dist < 1.0);

	rpo = m_angio.Position(seg.tip(1).pt);
	rp = seg.tip(1).pt.r;
	dist = distance(rpo.x, rp.x, rpo.y, rp.y, rpo.z, rp.z);
	assert(dist < 1.0);

	seg.m_nid = m_nsegs;
	SimulationTime& sim_time = m_angio.CurrentSimTime();
	seg.SetTimeOfBirth(sim_time.t);
	m_nsegs++;
	m_frag.push_front(seg);
	assert(m_nsegs == (int)m_frag.size());
}

//-----------------------------------------------------------------------------
// Determine the orientation of a newly created segment
vec3d Culture::FindGrowDirection(Segment::TIP& tip) const
{
    // Find the component of the new vessel direction determined by collagen fiber orientation    
    vec3d coll_dir = m_angio.CollagenDirection(tip.pt);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = (coll_dir*m_cultParams->vessel_orient_weights.x + per_dir*m_cultParams->vessel_orient_weights.y)/(m_cultParams->vessel_orient_weights.x+m_cultParams->vessel_orient_weights.y);
	new_dir.unit();

	return new_dir;
}

//-----------------------------------------------------------------------------
// Find the density-based length scale factor at a point of the grid
double Culture::FindDensityScale (const GridPoint& pt) const
{
	double coll_den = m_angio.FindECMDensity(pt);

	// Determine the density scaling factor using the function defined by a, b, c
	double den_scale;
	den_scale = m_cultParams->m_density_scale_factor.x + m_cultParams->m_density_scale_factor.y
		*exp(-m_cultParams->m_density_scale_factor.z*coll_den);

	if (den_scale < 0.0) den_scale = 0.0;

	return den_scale;
}

//-----------------------------------------------------------------------------
// creates a new segment to connect close segments
Segment Culture::ConnectSegment(Segment& it1, Segment& it2, int k1, int k2)
{
 	Segment seg;
 	seg.seed(it1.seed());
	seg.vessel(it1.vessel());
 	
	seg.tip(0).pt = it1.tip(k1).pt;
	seg.tip(1).pt = it2.tip(k2).pt;
	seg.Update();
 	
 	seg.tip(0).bactive = false;
 	seg.tip(1).bactive = false;

	return seg;
 }

//-----------------------------------------------------------------------------
// Updates the list of active tips.
void Culture::FindActiveTips()
{
	m_active_tips.clear();
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		if (it->tip(0).bactive) m_active_tips.push_back(&it->tip(0));
		if (it->tip(1).bactive) m_active_tips.push_back(&it->tip(1));
	} 
}

//-----------------------------------------------------------------------------
// Use the displacement field from the FE solution to update microvessels into the current configuration
void Culture::Update()
{
	// loop over all fragments
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		// Iterate through both segment tips
		for (int k=0; k<2; ++k)
		{
			// get the tip
			Segment::TIP& tip = it->tip(k);

			// Update position
			assert(tip.pt.nelem >= 0);
			tip.pt.r = m_angio.Position(tip.pt);
		}
		
		// Recalculate the segment's length and unit vector based on it's new position
		it->Update();
	}

	// Update the total vascular length within the simulation   
    m_total_length = 0.;
    for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
    {
        m_total_length += it->length();
    }
}   

void DirectionalWeights(double da, double dw[2])
{
	double degree_isotropy = 1.0-da;
	dw[0] = degree_isotropy;
	dw[1] = da;
}