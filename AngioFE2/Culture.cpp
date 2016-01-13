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
	m_anast_dist = 75.0;

	// parameter for growth curve (TODO: Maybe move all this to Culture?)
    m_a = 1900.0;
    m_b = 1.4549;
    m_x0 = 4.9474;
    m_y0 = -19.1278;
    m_d = m_y0 + m_a/(1.0 + exp(m_x0/m_b)); // Initial value of growth curve (t = 0)

	// Probability that initial fragments form branches
	m_init_branch_prob = 0.2697;

	m_vess_length = m_d;

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

	return true;
}

//-----------------------------------------------------------------------------
// Create initial fragments
void Culture::SeedFragments(SimulationTime& time)
{
	for (int i=0; i < m_ninit_frags; ++i)
	{
		// Create an initial segment
		Segment seg = createInitFrag();

		// Give the segment the appropriate label
		seg.m_nseed = i;

		// Set the segment vessel as the segment label
		seg.m_nvessel = seg.m_nseed;

		// add it to the list
		AddSegment(seg);
	}

	// set the number of vessels
    m_num_vessel = m_ninit_frags;

	// Update the active tip container
	FindActiveTips();
}

//-----------------------------------------------------------------------------
// Generates an initial fragment that lies inside the grid.
Segment Culture::createInitFrag()
{
	Grid& grid = m_angio.GetGrid();

	// Create a segment
    Segment seg;

	// Set seg length to value of growth function at t = 0
	double seg_length = m_vess_length;
    
	// We create only segments that lie inside the grid.
	// Since the creation of such a segment may fail (i.e. too close to boundary)
	// we loop until we find one.
	// TODO: This could potentially create an infinite loop so I should probably safeguard against it.
	GridPoint p1;
	do
	{
		// pick a random element
		// The factor 0.999 is because there is a small chance that
		// frand() returns 1 in which case elem_num would be invalid.
		int elem_num = int(0.999*frand()*grid.Elems());

		// generate random natural coordinates
		vec3d q = vrand();

		// set the position of the first tip
		seg.tip(0).pt = GridPoint(elem_num, q);
		seg.tip(0).rt = grid.Position(seg.tip(0).pt);
    
		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = CollagenDirection(seg.tip(0).pt);

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = seg.tip(0).rt + seg_vec*seg_length;

		// find the element where the second tip is
		grid.FindGridPoint(r1, p1);
	}
	while (p1.nelem == -1);

	// assign the grid point to the end tip
	seg.tip(1).pt = p1;
	seg.tip(1).rt = grid.Position(p1);
	
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
	return seg;
}

//-----------------------------------------------------------------------------
// Vessel elongation is represented by the addition of new line segments at the locations of the active sprouts.
void Culture::Grow(SimulationTime& time)
{
	// Determine length of new segments for this growth step        		
	UpdateNewVesselLength(time);

	// Elongate the active vessel tips
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		for (int k=0; k<2; ++k)
		{
			// only grow active tips
			if (it->tip(k).bactive)
			{
				// Create new vessel segment at the current tip existing segment
				Segment seg = GrowSegment(*it, k);

				// Add the segment to the network
				AddNewSegment(seg, k);
			}
	    }
    }

	if (yes_branching) Branch(time);

	// If anatomosis is turned on determine which segments form anatomoses
	if (yes_anast) Fuse(time);										

	// Remove buggy segments
	kill_dead_segs();

	// Locate all active growth tips
	FindActiveTips();
}

//-----------------------------------------------------------------------------
void Culture::SubGrowth(double scale)
{
	// Iterator through the list of active segment tips
	for (list<SegIter>::iterator tip_it = m_active_tips.begin(); tip_it != m_active_tips.end(); ++tip_it)
	{
		// Dereference the tip iterator to obtain the active segment
		Segment& seg = (*(*tip_it));

		// Step growth for the active segments
		if (seg.tip(0).bactive) seg.tip(0).rt = seg.tip(1).rt - seg.uvect()*(scale*seg.length());
		if (seg.tip(1).bactive) seg.tip(1).rt = seg.tip(0).rt + seg.uvect()*(scale*seg.length());
	}
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
// Create a new segment at the tip of an existing segment
Segment Culture::GrowSegment(Segment& it, int k, bool branch, bool bnew_vessel)
{
	// Make sure the tip is active
	assert(it.tip(k).bactive);

	// calculate the length scale factor based on collagen density
	double den_scale = FindDensityScale(it.tip(k).pt);

	// this is the new segment length
	double seg_length = den_scale*m_vess_length;

	// determine the growth direction
	vec3d seg_vec = FindDirection(it, it.tip(k).pt);

	// If new segment is a branch we modify the grow direction a bit
	if (branch)
	{
		vec3d coll_fib = CollagenDirection(it.tip(k).pt);
		seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
		seg_vec.unit();
	}

	// Create a new segment
	Segment seg;

	// transer seed label
	seg.m_nseed  = it.m_nseed;

	// assign vessel ID
	if (bnew_vessel)
	{
		seg.m_nvessel = m_num_vessel++;
		assert(seg.m_nvessel > it.m_nvessel);
	}
	else seg.m_nvessel = it.m_nvessel;

	// grow from the end point
	if (k == 1)
	{
		// the first point is a copy of the last point of the source segment
		seg.tip(0).rt = it.tip(1).rt;
		seg.tip(0).pt = it.tip(1).pt;

		// position the end point
		seg.tip(1).rt = seg.tip(0).rt + seg_vec*seg_length;

        seg.tip(1).bactive = true;		// Turn on end tip of new segment
		seg.tip(0).bactive = false;		// Turn off origin tip of new segment
	}
	else // grow from the start point
	{
		// the first point is a copy of the last point of the source segment
		seg.tip(1).rt = it.tip(0).rt;
		seg.tip(1).pt = it.tip(0).pt;
		
		// position the end point (notice negative sign)
		seg.tip(0).rt = seg.tip(1).rt - seg_vec*seg_length;

		seg.tip(0).bactive = true;		// Turn on end tip of new segment
		seg.tip(1).bactive = false;		// Turn off origin tip of new segment
	}

	// Turn off previous segment tip
	it.tip(k).bactive = false;

	// move the body force to the new tip
	seg.tip(k).bdyf_id = it.tip(k).bdyf_id;

	// update length and unit vector
	seg.Update();

	return seg;
}

//-----------------------------------------------------------------------------
// Branching phase of the growth step.
void Culture::Branch(SimulationTime& time)
{
	// Elongate the active vessel tips
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		// Make sure segment has not reached a boundary
		// and has not undegone anastimoses.
		if ((it->GetFlag(Segment::BC_DEAD)==false)&&
			(it->GetFlag(Segment::ANAST  )==false))
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

//			fprintf(stdout, "branching probability: %lg\n", bprob);

			// If branching is turned on determine if the segment forms a new branch
			if ((frand() < bprob) || it->GetFlag(Segment::INIT_BRANCH)) BranchSegment(*it, k);
		}
	}

	// Remove buggy segments
	kill_dead_segs();
}

//-----------------------------------------------------------------------------
// Branching is modeled as a random process, determine is the segment passed to this function forms a branch
void Culture::BranchSegment(Segment& it, int k)
{
	// Increase the number of branches.
	m_num_branches = m_num_branches + 1;

	// Turn off the initial branch flag
	it.SetFlagOff(Segment::INIT_BRANCH);
    			
	// we must reactive the tip
	// (it will be deactivated by CreateNewSeg)
	it.tip(k).bactive = true;

	// Create the new vessel segment (with branch flag and new_vessel flag to true)
	Segment seg = GrowSegment(it,k, true, true);

	// Create a new sprout force for the branch
	CreateBranchingForce(seg);

	// Add it to the culture
	AddNewSegment(seg, k);	                                                              
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
	
	if (seg.tip(0).bactive)
	{
		tip = seg.tip(0).rt;															// Obtain the position of the new tip
		sprout_vect = -seg.uvect();	// notice negative sign
		seg.tip(0).bdyf_id = m_angio.total_bdyf - 1;											// Assign the body force ID
	}
		
	if (seg.tip(1).bactive)
	{
		vec3d tip = seg.tip(1).rt;
		sprout_vect = seg.uvect();
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
// Fuse
void Culture::Fuse(SimulationTime& time)
{
    for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		check4anast(*it, 0, time);
        check4anast(*it, 1, time);
	} 
}

//-----------------------------------------------------------------------------
// check4anast
void Culture::check4anast(Segment& it, int k, SimulationTime& time)
{
    if (it.tip(k).bactive == false) return;
	
	if (it.GetFlag(Segment::ANAST)) return;
	
    for (SegIter it2 = m_frag.begin(); it2 != m_frag.end(); ++it2)
	{                                                           
	    double dist0 = (it.tip(k).rt - it2->tip(0).rt).norm(); dist0 *= dist0;
		double dist1 = (it.tip(k).rt - it2->tip(1).rt).norm(); dist1 *= dist1;
        anastomose(dist0, dist1, k, it, *it2, time);
    } 
}

//-----------------------------------------------------------------------------
// anastomose
void Culture::anastomose(double dist0, double dist1, int k, Segment& it1, Segment& it2, SimulationTime& time)
{
	// make sure neither segments are anastomosed.
	if ((it1.GetFlag(Segment::ANAST)) || (it2.GetFlag(Segment::ANAST))) return;
    
	// make sure neither segments sprout from the same initial fragment
	// TODO: why is this not allowed?
    if (it1.m_nseed == it2.m_nseed) return;

    int kk = -1;
    if (dist0 <= m_anast_dist)
        kk = 0;
    else if (dist1 <= m_anast_dist)
        kk = 1;
  
    if (kk == -1) return;
                                                
	Segment seg = ConnectSegment(it1,it2,k,kk, time);      // CULTURE.connectSegment(segment 1, segment 2, tip 1, tip 2, grid, data, frag list container)
					                                            // This function will create a segment between to two segments to complete the anastomosis
	AddSegment(seg);                                      // Append new segment at the top of the list 'frag'                  
	it1.tip(k ).bactive = false;		// Deactivate tip of segment 1 after anastomosis
	it2.tip(kk).bactive = false;		// Deactivate tip of segment 2 after anastomosis (tip-tip anastomosis only)
}

//-----------------------------------------------------------------------------
// Add a new segment to the culture.
// This will apply BCs to the new segment and may result in 
// the addition of several new segments. 
// k is the new tip
// (it is assumed that the other tip was form at the end of 
// another segment and is valid)
void Culture::AddNewSegment(Segment& seg, int k)
{
	// get the new tip
	Segment::TIP& new_tip = seg.tip(k);
	assert(new_tip.pt.nelem == -1);

	// Find the position of the new end point
	Grid& grid = m_angio.GetGrid();
	if (grid.FindGridPoint(new_tip.rt, new_tip.pt) == false)
	{
		// If we get here, the new end point lies outside the grid
		// In that case, we apply boundary conditions
		// (Note that this may create additional segments)
		bc.checkBC(seg, k);
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

	seg.m_nid = m_nsegs;
	seg.SetTimeOfBirth(m_angio.m_time.t);
	m_nsegs++;
	m_frag.push_front(seg);
	assert(m_nsegs == (int)m_frag.size());
}

//-----------------------------------------------------------------------------
// Determine the orientation of a newly created segment
vec3d Culture::FindDirection(Segment& it, GridPoint& pt)
{
    // Find the component of the new vessel direction determined by collagen fiber orientation    
    vec3d coll_dir = CollagenDirection(pt);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = it.uvect();

	vec3d new_dir = (coll_dir*m_W[0] + per_dir*m_W[3])/(m_W[0]+m_W[3]);
	new_dir.unit();

	return new_dir;
}

//-----------------------------------------------------------------------------
// Calculates the unit direction vector of the collagen fibers at a grid point.
vec3d Culture::CollagenDirection(GridPoint& pt)
{   
	// get the element
	Grid& grid = m_angio.GetGrid();
	Elem& elem = grid.ebin[pt.nelem];
        
    // Obtain shape function weights
    double shapeF[8];
    grid.shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);
    
    // Determine component of new vessel direction due to nodal collagen fiber orientation
	vec3d coll_angle = 
		((*elem.n1).m_collfib)*shapeF[0] + 
		((*elem.n2).m_collfib)*shapeF[1] + 
		((*elem.n3).m_collfib)*shapeF[2] + 
		((*elem.n4).m_collfib)*shapeF[3] + 
		((*elem.n5).m_collfib)*shapeF[4] + 
		((*elem.n6).m_collfib)*shapeF[5] + 
		((*elem.n7).m_collfib)*shapeF[6] + 
		((*elem.n8).m_collfib)*shapeF[7];

	// make unit vector
	coll_angle.unit();

    return coll_angle;
}

//-----------------------------------------------------------------------------
// Calculates the density scale factor at a grid point.
double Culture::FindDensityScale(const GridPoint& pt)
{
	Grid& grid = m_angio.GetGrid();

	// get the element
	Elem& elem = grid.ebin[pt.nelem];
        
    // Obtain shape function weights
    double shapeF[8];
    grid.shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);

	// evaluate the collagen density
	double coll_den = 
		shapeF[0]*(*elem.n1).m_ecm_den + 
		shapeF[1]*(*elem.n2).m_ecm_den + 
		shapeF[2]*(*elem.n3).m_ecm_den + 
		shapeF[3]*(*elem.n4).m_ecm_den + 
		shapeF[4]*(*elem.n5).m_ecm_den + 
		shapeF[5]*(*elem.n6).m_ecm_den + 
		shapeF[6]*(*elem.n7).m_ecm_den + 
		shapeF[7]*(*elem.n8).m_ecm_den;
    
	// scale the density
	double den_scale = ScaleDensity(coll_den);
	
	// all done
	return den_scale;
}

//-----------------------------------------------------------------------------
// TODO: Make a,b,c user settings.
double Culture::ScaleDensity(double coll_den)
{
	if (coll_den == 3.0) return 1.0;

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
Segment Culture::ConnectSegment(Segment& it1, Segment& it2, int k, int kk, SimulationTime& time)
{
 	Segment seg;
 	seg.m_nseed   = it1.m_nseed;
 	seg.m_nvessel = it1.m_nvessel;
 	
 	seg.tip(0).rt = it1.tip(k).rt;
	seg.tip(0).pt = it1.tip(k).pt;
	seg.tip(1).rt = it2.tip(kk).rt;
	seg.tip(1).pt = it2.tip(kk).pt;
	seg.Update();
 	

 	seg.tip(0).bactive = false;
 	seg.tip(1).bactive = false;
	seg.SetFlagOn(Segment::ANAST);
	it1.SetFlagOn(Segment::ANAST);
	it2.SetFlagOn(Segment::ANAST);

 	++m_num_anastom;
 	
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


double findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3])
{
	double t1[3], t2[3];
	double A, B, C, D, E, F, min_dist;
	double lambda, mu;
	
	t1[0] = b[0] - a[0];
	t1[1] = b[1] - a[1];
	t1[2] = b[2] - a[2];

	t2[0] = d[0] - c[0];
	t2[1] = d[1] - c[1];
	t2[2] = d[2] - c[2];

	A = vec_dot(t1,t1);
	B = 2*(vec_dot(a,t1) - vec_dot(t1,c));
	C = 2*(vec_dot(t1,t2));
	D = 2*(vec_dot(t2,c) - vec_dot(t2,a));
	E = vec_dot(t2,t2);
	F = vec_dot(a,a)+ vec_dot(c,c);

	if (4*A*E-C*C != 0)
	{
		min_dist = (B*C*D + B*B*E + C*C*F + A*(D*D-4*E*F))/(C*C - 4*A*E);
		if (min_dist <= 7)
		{
			mu = (2*A*D + B*C)/(C*C - 4*A*E);
			lambda = (C*mu - B)/(2*A);
			intersectpt[0] = a[0] + lambda*t1[0];
			intersectpt[1] = a[1] + lambda*t1[1];
			intersectpt[2] = a[2] + lambda*t1[2];
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
// Description: checks for intersection between a passed segment and all other existing segments that are not members
// of vessel containing the segment
void Culture::CheckForIntersection(Segment &seg, Segment& it)
{
	double p1[3], p2[3], pp1[3], pp2[3]; //tip points for segment to check intersection
	list<Segment>::iterator it2; //loop over segments in list
	double intersectpt[3]; //the intersection pt of vectors
	
	intersectpt[0] = 0;
	intersectpt[1] = 0;
	
	double lambda;

	vec3d& r0 = seg.tip(0).rt;
	vec3d& r1 = seg.tip(1).rt;
	
	if (seg.tip(1).bactive)
	{
		p1[0] = r0.x;
		p2[0] = r1.x;
		p1[1] = r0.y;
		p2[1] = r1.y;
		p1[2] = r0.z;
		p2[2] = r1.z;
	}
	else if (seg.tip(0).bactive)
	{
		p1[0] = r1.x;
		p2[0] = r0.x;
		p1[1] = r1.y;
		p2[1] = r0.y;
		p1[2] = r1.z;
		p2[2] = r0.z;
	}

	Culture& cult = m_angio.GetCulture();
	for (it2 = cult.m_frag.begin(); it2 != cult.m_frag.end(); ++it2)
	{
		if (it.m_nseed!=it2->m_nseed)
		{
			pp1[0] = r0.x;
			pp1[1] = r0.y;
			pp1[2] = r0.z;
			pp2[0] = r1.x;
			pp2[1] = r1.y;
			pp2[2] = r1.z;

			lambda = findIntersect(p1,p2,pp1,pp2,intersectpt);
			if (lambda >=0 && lambda <=1)
			{
				if (seg.tip(1).bactive)
				{
					seg.tip(1).rt.x = intersectpt[0];
					seg.tip(1).rt.y = intersectpt[1];
					seg.tip(1).rt.z = intersectpt[2];
					seg.tip(1).bactive = false;
				}
				else if (seg.tip(0).bactive)
				{
					seg.tip(0).rt.x = intersectpt[0];
					seg.tip(0).rt.y = intersectpt[1];
					seg.tip(0).rt.z = intersectpt[2];
					seg.tip(0).bactive = false;
				}
				cout << "3D intersection" << endl;
				++m_num_anastom;
			}
		}
	}


	return;
}

//-----------------------------------------------------------------------------
// Remove dead segments.
// TODO: It think this is a sloppy way of sweeping bugs under the carpet. Remove this.
void Culture::kill_dead_segs()
{  
	if (m_angio.kill_off == false){
		list<Segment>::iterator it;
    
		for (it = m_frag.begin(); it != m_frag.end(); ++it){
			if (it->mark_of_death == true){
				vec3d& r0 = it->tip(0).rt;
				vec3d& r1 = it->tip(1).rt;
				it = m_frag.erase(it);
				assert(false);
			}
		}
	}
}


//-----------------------------------------------------------------------------
// Updates the list of active tips.
// TODO: I think this logic is flawed. This presumes that only one tip is active
//       yet the initial fragments have both tips active.
void Culture::FindActiveTips()
{
	list<Segment>::iterator it;
           
	m_active_tips.clear();
	
	for (it = m_frag.begin(); it != m_frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		if (it->tip(0).bactive)
			m_active_tips.push_back(it);
		else if (it->tip(1).bactive)
			m_active_tips.push_back(it);
	} 
			
    return;
}
