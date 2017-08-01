#include "StdAfx.h"
#include "Culture.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"
#include "GrowDirectionModifier.h"
#include "FECore/FEMesh.h"
#include <random>
#include <regex>
#include <map>
#define _USE_MATH_DEFINES

void DirectionalWeights(double da, double dw[2]);

std::vector<double> SegTo1TipPos(Segment * seg)
{
	std::vector<double> rv;

	const vec3d & pos = seg->tip(1).pos();

	rv.emplace_back(pos.x);
	rv.emplace_back(pos.y);
	rv.emplace_back(pos.z);

	return rv;
}
std::vector<double> units(3, 1.0);

//-----------------------------------------------------------------------------
Culture::Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp, FragmentBranching *fbr) : m_angio(angio), tips(SegTo1TipPos, ndim_distance, ndim_distance_to_plane, units)
{
	assert(matl && cp);
	m_cultParams = cp;

	m_num_zdead = 0;

	// TODO: I think it would make sense to tie the anastomosis distance 
	//       to the vessel diameter.
    m_total_length = 0.;
	m_num_vessel = 0;
	m_num_branches = 0;		// Initialize branching counter

	// Initialize counters
	m_num_anastom = 0;
	m_nsegs = 0;			// Initialize segment counter
	m_pmat = matl;
	fbrancher = fbr;
}

//-----------------------------------------------------------------------------
Culture::~Culture()
{
}

//-----------------------------------------------------------------------------
// Initialize the culture
bool Culture::Init()
{
	if (!m_pmat->bc)
		return false;
	m_pmat->bc->SetCulture(this);

	//intialize the Fragment Seeder
	m_pmat->fseeder->SetCulture(this);

	fbrancher->SetCulture(this);

	m_pmat->gdms->SetCulture(this);

	// do the initial seeding
	if (!m_pmat->fseeder->SeedFragments(m_angio.CurrentSimTime(), this))
		return false;

	return true;
}

//this formula is from the previous version of the plugin
//this is derived from the sigmoid curve http://www.sci.utah.edu/publications/edgar13/Edgar_CMBBE2013.pdf
double Culture::SegmentLength(double starttime, double grow_time) const
{
	//currently this ignores the starttime this aligns with the other constraints on step refinement
	double lc = m_cultParams->growth_length_over_time * grow_time;

	return lc*m_cultParams->m_length_adjustment;
}

//-----------------------------------------------------------------------------
// Create a new segment at the (active) tip of an existing segment
Segment Culture::GrowSegment(Segment::TIP& tip, double start_time, double grow_time, bool branch, bool bnew_vessel)
{

	// Make sure the tip is active
	assert(tip.bactive);

	// this is the new segment length
	double seg_length = 0.0;
	
	// determine the growth direction
	vec3d seg_vec = m_pmat->gdms->ApplyModifiers( vec3d(), tip, m_pmat, branch, start_time, grow_time, seg_length);


	// Create a new segment
	Segment seg;

	seg.expected_length = seg_length;

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

	seg.tip(0).connected = tip.parent;

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
// Branching is modeled as a random process, determine is the segment passed to this function forms a branch
void Culture::BranchSegment(Segment::TIP& tip, double starttime, double grow_time)
{
	// we must reactive the tip
	// (it will be deactivated by GrowSegment)
	tip.bactive = true;

	//calculate the density gradinet if above the threshold set the grow direction
	
	m_num_branches++;
	

	Segment seg = GrowSegment(tip, starttime, grow_time ,true, true);

	//copy the branch distance in from the tip
	seg.tip(0).length_to_branch = tip.length_to_branch;
	seg.tip(1).wait_time_to_branch = tip.wait_time_to_branch;
	seg.tip(0).wait_time_to_branch = tip.wait_time_to_branch;

	//will not scale with multiple threads
	double prev_min_seg_length = m_cultParams->min_segment_length;
	m_cultParams->min_segment_length = 0.0;
	if(seg.length() != 0.0)
	{
		// Add it to the culture
		AddNewSegment(seg);
	}
	else
	{
		ClearRecents();
	}
	

	//restore old min seg length
	m_cultParams->min_segment_length = prev_min_seg_length;
}

//-----------------------------------------------------------------------------
// Create a new sprout force component for a newly formed branch
// TODO: What if both tips are active?
void Culture::CreateBranchingForce(Segment& seg)
{															// Obtain the total number of sprouts
	m_angio.total_bdyf = m_pmat->Sprouts();

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

	m_pmat->AddSprout(tip, sprout_vect, stip.pt.ndomain, stip.pt.elemindex);
	m_angio.total_bdyf = m_pmat->Sprouts();
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
	recents.clear();
	//adding zero length segments shoudl be avoided
	//this will clear and refill recents
	assert(seg.length() > 0.0);
	assert(seg.tip_c(0).pt.nelem >= 0);

	// get the new tip
	Segment::TIP& new_tip = seg.tip(1);
	//assert(new_tip.pt.nelem == -1);
	//init done elsewhere
	assert(new_tip.bactive);

	m_pmat->bc->CheckBC(seg);
}

void Culture::AddNewSegmentNoClear(Segment& seg)
{
	//adding zero length segments shoudl be avoided
	//this will clear and refill recents
	assert(seg.length() > 0.0);
	assert(seg.tip_c(0).pt.nelem >= 0);
	assert(seg.tip_c(0).connected);


	// get the new tip
	Segment::TIP& new_tip = seg.tip(1);
	//assert(new_tip.pt.nelem == -1);
	//init done elsewhere
	assert(new_tip.bactive);

	m_pmat->bc->CheckBC(seg);
}

//-----------------------------------------------------------------------------
// This function adds a segment to the list.
// It also updates the segment counter and assigns Segment::m_nid and sets the time of birth.
// Note that we also push to the front of the list so that we don't corrupt
// any active iterators.
void Culture::AddSegment(Segment& seg)
{
	// let's check a few things
	assert(seg.tip_c(0).pt.nelem >= 0);
	assert(seg.tip_c(1).pt.nelem >= 0);

	assert(seg.tip_c(0).pt.ndomain != nullptr);
	assert(seg.tip_c(1).pt.ndomain != nullptr);

	assert(seg.tip_c(0).pt.elemindex >= 0);
	assert(seg.tip_c(1).pt.elemindex >= 0);

	assert(seg.tip_c(0).nseed >= 0);
	assert(seg.tip_c(1).nseed >= 0);

	assert(seg.tip_c(0).nvessel >= 0);
	assert(seg.tip_c(1).nvessel >= 0);

	assert(seg.seed() >= 0);
	assert(seg.vessel() >= 0);

	assert(seg.length() > 0.0);

	assert(seg.length() < seg.expected_length*1.01);

	vec3d rpo = m_angio.Position(seg.tip(0).pt);
	vec3d rp = seg.tip(0).pt.r;
	double dist = distance(rpo.x, rp.x, rpo.y, rp.y, rpo.z, rp.z);
	assert(dist < 1.0);

	rpo = m_angio.Position(seg.tip(1).pt);
	rp = seg.tip(1).pt.r;
	dist = distance(rpo.x, rp.x, rpo.y, rp.y, rpo.z, rp.z);
	assert(dist < 1.0);

	seg.m_nid = m_nsegs;
	m_nsegs++;
	m_frag.emplace_front(seg);
	assert(m_nsegs == (int)m_frag.size());
	
	//add the segment that was just created
	recents.emplace_back(&m_frag.front());
	segments_per_step.emplace_back(&m_frag.front());

	tips.insert(&m_frag.front());
	//add the sprouts to the material
	if (seg.tip(0).bactive)
	{
		m_pmat->AddSprout(seg.tip(0));
	}
	if (seg.tip(1).bactive)
	{
		m_pmat->AddSprout(seg.tip(1));
	}
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
		if (it->tip(0).bactive) m_active_tips.emplace_back(&it->tip(0));
		if (it->tip(1).bactive) m_active_tips.emplace_back(&it->tip(1));
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
	//rebuild the kd tree as the positions withn it may have moved
	tips.rebuild();
}   