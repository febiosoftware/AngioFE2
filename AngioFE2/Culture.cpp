#include "stdafx.h"
#include <iostream>
#include "Culture.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"

//-----------------------------------------------------------------------------
Culture::Culture(FEAngio& angio) : bc(angio), m_angio(angio)
{
	m_W[0] = 0;
	m_W[1] = 0;
	m_W[2] = 0;
	m_W[3] = 0;

	m_ninit_frags = 0;
    m_num_vessel = 0;
	m_num_branches = 0;		// Initialize branching counter

	m_num_zdead = 0;

	// TODO: I think it would make sense to tie the anastomosis distance 
	//       to the vessel diameter.
	m_anast_dist = 8.6;

    m_total_length = 0.;
    
	// parameter for growth curve
    m_a = 1900.0;
    m_b = 1.4549;
    m_x0 = 4.9474;
    m_y0 = -19.1278;
    double d = m_y0 + m_a/(1.0 + exp(m_x0/m_b)); // Initial value of growth curve (t = 0)

	m_vess_length = d;

	m_init_length = 0.0;

	// Probability that initial fragments form branches
	m_init_branch_prob = 0.2697;

	m_length_adjust = 1.0;

	// Flags
	yes_branching = true;
	yes_anast = true;

	// branch probability
	m_branch_chance = 0.1;

	// Initialize counters
	m_num_anastom = 0;
	m_nsegs = 0;			// Initialize segment counter
}

//-----------------------------------------------------------------------------
Culture::~Culture()
{

}

//-----------------------------------------------------------------------------
// Initialize the culture
bool Culture::Init()
{
	// If branching is turned off, we set the branching chance to zero
	if (yes_branching == false) m_branch_chance = 0.0;

	// make sure the initial length is initialized 
	if (m_init_length <= 0.0) m_init_length = m_vess_length;

	// do the initial seeding
	return SeedFragments(m_angio.CurrentSimTime());
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool Culture::SeedFragments(SimulationTime& time)
{
	for (int i=0; i < m_ninit_frags; ++i)
	{
		// Create an initial segment
		Segment seg;
		if (createInitFrag(seg) == false) return false;

		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		// add it to the list
		AddSegment(seg);
	}

	// init vessel counter
    m_num_vessel = m_ninit_frags;

	// Update the active tip container
	FindActiveTips();

	return true;
}

//-----------------------------------------------------------------------------
// Generates an initial fragment that lies inside the grid.
bool Culture::createInitFrag(Segment& seg)
{
	Grid& grid = m_angio.GetGrid();

	// Set seg length to value of growth function at t = 0
	double seg_length = m_init_length;
    
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
		int elem_num = int(0.999*frand()*grid.Elems());

		// generate random natural coordinates
		vec3d q = vrand();

		// set the position of the first tip
		p0 = grid.FindGridPoint(elem_num, q);
    
		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = grid.CollagenDirection(p0);

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		grid.FindGridPoint(r1, p1);

		ntries++;
	}
	while ((p1.nelem == -1) && (ntries < MAX_TRIES));

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
	if (frand() < m_init_branch_prob) seg.SetFlagOn(Segment::INIT_BRANCH);
	
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
	if (yes_branching) BranchVessels(time);

	// Anastomosis phase
	if (yes_anast) FuseVessels();

	// Update all active growth tips
	FindActiveTips();
}

//-----------------------------------------------------------------------------
// The new vessel lenght (i.e. the lenght of new vessels) is a function
// of time. This function is called at the start of each growth step (Culture::Grow)
// and updates the new vessel length.
void Culture::UpdateNewVesselLength(SimulationTime& time)
{
	double lc = m_a/(1.+ exp(-(time.t - m_x0)/m_b));
	lc -= m_a/(1. + exp(-(time.t - time.dt-m_x0)/m_b));
	
	m_vess_length = lc*m_length_adjust;
}

//-----------------------------------------------------------------------------
// Grow phase
void Culture::GrowVessels()
{
	for (TipIter it = m_active_tips.begin(); it != m_active_tips.end(); ++it)
	{
		Segment::TIP& tip = *(*it);

		// Create new vessel segment at the current tip existing segment
		Segment seg = GrowSegment(tip);

		// Add the segment to the network
		// This will also enforce the boundary conditions
		AddNewSegment(seg);
    }
}

//-----------------------------------------------------------------------------
// Create a new segment at the (active) tip of an existing segment
Segment Culture::GrowSegment(Segment::TIP& tip, bool branch, bool bnew_vessel)
{
	Grid& grid = m_angio.GetGrid();

	// Make sure the tip is active
	assert(tip.bactive);

	// calculate the length scale factor based on collagen density
	double den_scale = FindDensityScale(tip.pt);

	// this is the new segment length
	double seg_length = den_scale*m_vess_length;

	// determine the growth direction
	vec3d seg_vec = FindGrowDirection(tip);

	// If new segment is a branch we modify the grow direction a bit
	if (branch)
	{
		// TODO: what's the logic here? Why the 0.5 factor?
		//      If the vessel is aligned with the collagen (and the initial fragments are)
		//      then  the new branch will overlap the old segment. 
		vec3d coll_fib = grid.CollagenDirection(tip.pt);
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

	// move the body force to the new tip
	seg.tip(1).bdyf_id = tip.bdyf_id;

	// Turn off previous segment tip
	tip.bactive = false;

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
			double bprob = den_scale*dt*m_branch_chance/t;
			assert(bprob < 1.0);

			// If branching is turned on determine if the segment forms a new branch
			if (frand() < bprob)
			{
				// Turn off the initial branch flag
				it->SetFlagOff(Segment::INIT_BRANCH);
				BranchSegment(it->tip(k));

				// Increase the number of branches.
				m_num_branches = m_num_branches + 1;
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

	// Create a new sprout force for the branch
	CreateBranchingForce(seg);

	// Add it to the culture
	AddNewSegment(seg);
}

//-----------------------------------------------------------------------------
// Create a new sprout force component for a newly formed branch
// TODO: What if both tips are active?
void Culture::CreateBranchingForce(Segment& seg)
{
	m_angio.total_bdyf = 0;																// Obtain the total number of sprouts
	if (m_angio.m_pmat) m_angio.total_bdyf = m_angio.m_pmat->Sprouts();
	else m_angio.total_bdyf = m_angio.m_pbf->Sprouts();

	vec3d tip, sprout_vect;
	vec3d seg_vec = seg.uvect();
	
	if (seg.tip(0).bactive)
	{
		tip = seg.tip(0).pos();															// Obtain the position of the new tip
		sprout_vect = -seg_vec;	// notice negative sign
		seg.tip(0).bdyf_id = m_angio.total_bdyf - 1;											// Assign the body force ID
	}
		
	if (seg.tip(1).bactive)
	{
		vec3d tip = seg.tip(1).pos();
		sprout_vect = seg_vec;
		seg.tip(1).bdyf_id = m_angio.total_bdyf - 1;
	}

	if (m_angio.m_pmat)
	{
		m_angio.m_pmat->AddSprout(tip, sprout_vect);
		m_angio.total_bdyf = m_angio.m_pmat->Sprouts();
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
							if (dist < m_anast_dist)
							{
								// create a segment between the two segments to complete the anastomosis
								Segment seg = ConnectSegment(*it1, *it2, k1, k2);

								// mark segments
								seg.SetFlagOn(Segment::ANAST);
								it1->SetFlagOn(Segment::ANAST);
								it2->SetFlagOn(Segment::ANAST);

								// add it to the list
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
	assert(seg.tip(0).pt.nelem >= 0);

	// get the new tip
	Segment::TIP& new_tip = seg.tip(1);
	assert(new_tip.pt.nelem == -1);
	assert(new_tip.bactive);

	// Find the position of the new end point
	Grid& grid = m_angio.GetGrid();
	if (grid.FindGridPoint(new_tip.pos(), new_tip.pt) == false)
	{
		// If we get here, the new end point lies outside the grid
		// In that case, we apply boundary conditions
		// (Note that this may create additional segments)
		bc.CheckBC(seg);
	}
	else
	{
		// everything looks good, so let's add the segment
		AddSegment(seg);
	}
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

	assert(seg.tip(0).nseed >= 0);
	assert(seg.tip(1).nseed >= 0);

	assert(seg.tip(0).nvessel >= 0);
	assert(seg.tip(1).nvessel >= 0);

	assert(seg.seed() >= 0);
	assert(seg.vessel() >= 0);

	assert(seg.length() < 1.1*m_vess_length);

	seg.m_nid = m_nsegs;
	SimulationTime& sim_time = m_angio.CurrentSimTime();
	seg.SetTimeOfBirth(sim_time.t);
	m_nsegs++;
	m_frag.push_front(seg);
	assert(m_nsegs == (int)m_frag.size());
}

//-----------------------------------------------------------------------------
// Determine the orientation of a newly created segment
vec3d Culture::FindGrowDirection(Segment::TIP& tip)
{
	Grid& grid = m_angio.GetGrid();

    // Find the component of the new vessel direction determined by collagen fiber orientation    
    vec3d coll_dir = grid.CollagenDirection(tip.pt);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = (coll_dir*m_W[0] + per_dir*m_W[3])/(m_W[0]+m_W[3]);
	new_dir.unit();

	return new_dir;
}

//-----------------------------------------------------------------------------
// Find the density-based length scale factor at a point of the grid
// TODO: Make a,b,c user settings.
double Culture::FindDensityScale(const GridPoint& pt)
{
	Grid& grid = m_angio.GetGrid();
	double coll_den = grid.FindECMDensity(pt);
    
	//if (coll_den == 3.0) return 1.0;

	// Determine the density scaling factor using the function defined by a, b, c
	double den_scale;
	double a = -0.016;
	double b = 5.1605;
	double c = 0.5112;
	den_scale = a + b*exp( -c*coll_den );

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
// Helper function for finding the intersection between two lines.
// variables:
// a : start x,y point on segment1
// b : end x,y point on segment1
// c : start x,y point on segment2
// d : end x,y point on segment2
// t1 : displacement vector for segment1 (b-a)
// t2 : displacement vector for segment2 (d-c)
// lambda : distance to move along t1 from its start point
// mu : distance to move along t2 from its start point

// General case for the minimum distance between two lines
// d^2 = [BCD + B^2E + C^2F + A(D^2-4EF)]/[C^2 - 4AE]
// where
// A = t1.t1
// B = 2(a.t1 - t1.c)
// C = 2(t1.t2)
// D = 2(t2.c - t2.a)
// E = t2.t2
// F = a.a + c.c
// This will have a unique solution if
// | 2A -C |
// | -C	2E | is non-zero
// closest point on t1 is : a + lambda*t1
// where
// lambda = (C*mu - B)/(2A)
// and
// mu = (2*A*D + B*C)/(C^2 - 4*A*E) 


double findIntersect(vec3d& a, vec3d& b, vec3d& c, vec3d& d, vec3d& intersectpt)
{
	double lambda;
	
	vec3d t1 = b - a;
	vec3d t2 = d - c;

	double A = t1* t1;
	double B = 2*(t1*a) - t1*c;
	double C = 2*(t1*t2);
	double D = 2*(t2*c - t2*a);
	double E = t2*t2;
	double F = a*a+ c*c;

	if (4*A*E-C*C != 0)
	{
		double min_dist = (B*C*D + B*B*E + C*C*F + A*(D*D-4*E*F))/(C*C - 4*A*E);
		if (min_dist <= 7)
		{
			double mu = (2*A*D + B*C)/(C*C - 4*A*E);
			lambda = (C*mu - B)/(2*A);
			intersectpt = a + t1*lambda;
		}
		else
		{
			return 1000;
		}
	}
	else
	{
		return 1000;
	}

	return lambda;
}

//-----------------------------------------------------------------------------
// Checks for intersection between a passed segment and all other existing segments that are not members
// of vessel containing the segment
// TODO: This is currently not used but might be useful for doing segment-segment intersections.
void Culture::CheckForIntersection(Segment &seg, Segment& it)
{
	vec3d p1, p2, pp1, pp2; //tip points for segment to check intersection
	vec3d intersectpt; //the intersection pt of vectors
	
	const vec3d& r0 = seg.tip(0).pos();
	const vec3d& r1 = seg.tip(1).pos();
	
	if (seg.tip(1).bactive)
	{
		p1 = r0;
		p2 = r1;
	}
	else if (seg.tip(0).bactive)
	{
		p1 = r1;
		p2 = r0;
	}

	for (SegIter it2 = m_frag.begin(); it2 != m_frag.end(); ++it2)
	{
		if (it.seed() !=it2->seed())
		{
			pp1 = r0;
			pp2 = r1;

			double lambda = findIntersect(p1,p2,pp1,pp2,intersectpt);
			if (lambda >=0 && lambda <=1)
			{
				if (seg.tip(1).bactive)
				{
					seg.tip(1).pt.r = intersectpt;
					seg.tip(1).bactive = false;
				}
				else if (seg.tip(0).bactive)
				{
					seg.tip(0).pt.r = intersectpt;
					seg.tip(0).bactive = false;
				}
				cout << "3D intersection" << endl;
				++m_num_anastom;
			}
		}
	}
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
	Grid& grid = m_angio.GetGrid();

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
			tip.pt.r = grid.Position(tip.pt);
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
