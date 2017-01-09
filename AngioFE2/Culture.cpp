#include "StdAfx.h"
#include "Culture.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include "angio3d.h"
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"
#include <random>
#include <regex>
#include <map>
#define _USE_MATH_DEFINES

void DirectionalWeights(double da, double dw[2]);

CultureParameters::CultureParameters()
{
	m_boundary_condition_type[0] = 's';
	m_boundary_condition_type[1] = '\0';
	vessel_file[0] = '\0';
}

//-----------------------------------------------------------------------------
Culture::Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp) : m_angio(angio)
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
	case 2:
		fseeder = new MDByVolumeFragmentSeeder(m_cultParams, m_angio);
		break;
	case 3:
		fseeder = new MDAngVessFileFragmentSeeder(m_cultParams, m_angio);
		break;
	default:
		assert(false);
	}

	switch (m_cultParams->branching_scheme)
	{
	case 0:
		fbrancher = new NoFragmentBranching(this);
		break;
		
	case 1:
		fbrancher = new ForwardFragmentBranching(this);
		break;
		/* removed until fully implemeented
	case 2:
		fbrancher = new PsuedoDeferedFragmentBranching(this);
		break;
		*/
	default:
		assert(false);
	}


	//intialize all of the grow direction modifiers
	if (m_cultParams->m_anastomosis)
	{
		grow_direction_modifiers.push_back(new AnastamosisGrowDirectionModifier(this, 1));
	}
	//this should always be the first thing in the queue
	grow_direction_modifiers.push_back(new DefaultGrowDirectionModifier(this, -1));

	//consider changing the priorities if needed / let the user set the priorities
	grow_direction_modifiers.push_back(new GradientGrowDirectionModifier(this, 0));

	grow_direction_modifiers.push_back(new BranchGrowDirectionModifier(this, 2));

	std::sort(grow_direction_modifiers.begin(), grow_direction_modifiers.end(), [](GrowDirectionModifier * a, GrowDirectionModifier * b)
	{
		return a->Priority() < b->Priority();
	});


	// do the initial seeding
	if (!fseeder->SeedFragments(m_angio.CurrentSimTime(), this))
		return false;

	return true;
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
// the classic mode still uses a biased growth direction
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
		sgi.domain = &mesh->Domain(culture->m_pmat->domains[ddist(m_angio.rengine)]);
		sgi.ielement = edist(m_angio.rengine);
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
bool FragmentSeeder::createInitFrag(Segment& seg, SegGenItem & item, Culture * culture)
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
		vec3d q = m_angio.uniformInUnitCube();

		// set the position of the first tip
		p0 = m_angio.FindGridPoint(item.domain, item.ielement, q);

		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = m_angio.uniformRandomDirection();
		seg_vec.unit();

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		culture->m_pmat->FindGridPoint(r1,p1);
		//need to check the domain is legal
		ntries++;
	} while (((p1.nelem == -1) || (std::find(domains.begin(), domains.end(), p1.ndomain) == domains.end())) && (ntries < MAX_TRIES));

	if (p1.nelem == -1)  return false;

	// assign the grid points
	seg.tip(0).pt = p0;
	seg.tip(1).pt = p1;

	// update length and unit vector
	seg.Update();

	// make both tips active
	seg.tip(0).bactive = true;
	seg.tip(1).bactive = true;

	seg.SetTimeOfBirth(0.0);
	seg.tip(0).length_to_branch = culture->fbrancher->GetLengthToBranch();
	seg.tip(1).length_to_branch = culture->fbrancher->GetLengthToBranch();

	assert(seg.tip(0).length_to_branch > 0.0);
	assert(seg.tip(1).length_to_branch > 0.0);

	// Mark segment as an initial fragment
	seg.SetFlagOn(Segment::INIT_SPROUT);

	// all good
	return true;
}

static size_t findElement(double val, int lo, int high, double * begin, double * end)
{
	int mid = lo + (high - lo) / 2;
	if (val < begin[mid])
	{
		return findElement(val, lo, mid -1, begin, end);
	}
	else if (val > end[mid])
	{
		return findElement(val, mid +1, high, begin, end);
	}
	else
	{
		return mid;
	}
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool MDByVolumeFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	FEMesh * mesh = m_angio.GetMesh();
	int elementsInDomains = 0;
	for (size_t i = 0; i < culture->m_pmat->domains.size(); i++)
	{
		elementsInDomains += m_angio.GetMesh()->Domain(culture->m_pmat->domains[i]).Elements();
		domains.push_back(&m_angio.GetMesh()->Domain(culture->m_pmat->domains[i]));
	}
	double * totalWeightsBegin = new double[elementsInDomains];
	double * totalWeightsEnd = new double[elementsInDomains];
	FEElement ** elements = new FEElement*[elementsInDomains];
	int k = 0;
	double cweight = 0.0;
	for (size_t i = 0; i < domains.size(); i++)
	{
		for (int j = 0; j < domains[i]->Elements(); j++)
		{
			FEElement & el = domains[i]->ElementRef(j);
			elements[k] = &el;
			totalWeightsBegin[k] = cweight;
			cweight += mesh->ElementVolume(el);
			totalWeightsEnd[k] = cweight;
			k++;
		}
	}
	std::uniform_real_distribution<double> voluchoice(0, cweight);
	SegGenItem sgi;

	for (int i = 0; i < culture_params->m_number_fragments; ++i)
	{
		//do a binary search to find the element that contains the volume
		double vol = voluchoice(m_angio.rengine);
		int ei = findElement(vol, 0, elementsInDomains - 1, totalWeightsBegin, totalWeightsEnd);
		// Create an initial segment
		Segment seg;
		FEElement * elem = elements[ei];
		sgi.ielement = elem->GetID() -1 - culture->m_pmat->meshOffsets[elem->GetDomain()];
		sgi.domain = elem->GetDomain();
		
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

	delete[] totalWeightsBegin;
	delete[] totalWeightsEnd;
	delete[] elements;
	return true;
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool MDAngVessFileFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	FEMesh * mesh = m_angio.GetMesh();
	infile.open(culture_params->vessel_file);
	if (infile.fail())
	{
		printf("error while opening input file\n");
		return false;
	}
	std::string line;
	int i = 1;
	vec3d p1, p2;
	std::regex whitespace("\\s+");
	while (std::getline(infile, line))
	{
		
		//any lines starting with a nondigit are considered comments
		if (!isdigit(line[0]))
		{
			continue;
		}
		int festate=0;
		double length = 0.0f, segtime = 0;
		line = regex_replace(line, whitespace, " ");
		if (9 != sscanf(line.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf", &festate, &segtime, &p1.x, &p1.y, &p1.z, &p2.x, &p2.y, &p2.z, &length))
		{
			//improperly formatted line
			//this line is skipped
			continue;
		}
			
		// Create an initial segment
		Segment seg;
		

		//find the elements that contain the point
		bool rv = culture->m_pmat->FindGridPoint(p1, seg.tip(0).pt);
		assert(rv);
		rv = culture->m_pmat->FindGridPoint(p2, seg.tip(1).pt);
		assert(rv);
		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		seg.Update();
		seg.tip(0).bactive = true;
		seg.tip(1).bactive = true;
		// add it to the list
		culture->AddSegment(seg);
		i++;
	}
	infile.close();

	// init vessel counter
	culture->m_num_vessel = culture_params->m_number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	return true;
}

//this formula is from the previous version of the plugin
double Culture::SegmentLength(double starttime, double grow_time) const
{
	double lc = m_cultParams->m_culture_a / (1. + exp(-(starttime - m_cultParams->m_x0) / m_cultParams->m_culture_b));
	lc -= m_cultParams->m_culture_a / (1. + exp(-(starttime - grow_time - m_cultParams->m_x0) / m_cultParams->m_culture_b));

	return lc*m_cultParams->m_length_adjustment;
}

//-----------------------------------------------------------------------------
// Create a new segment at the (active) tip of an existing segment
Segment Culture::GrowSegment(Segment::TIP& tip, double starttime, double grow_time, bool branch, bool bnew_vessel)
{

	// Make sure the tip is active
	assert(tip.bactive);

	// calculate the length scale factor based on collagen density
	double den_scale = FindDensityScale(tip.pt);

	// this is the new segment length
	double seg_length = SegmentLength(starttime, grow_time);

	// determine the growth direction
	vec3d seg_vec;
	//now run it through the different filters
	for (size_t i = 0; i < grow_direction_modifiers.size(); i++)
	{
		seg_vec = grow_direction_modifiers[i]->GrowModifyGrowDirection(seg_vec, tip, branch);
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
	seg.tip(1).length_to_branch = tip.length_to_branch;

	//will not scale with multiple threads
	double prev_min_seg_length = m_cultParams->min_segment_length;
	m_cultParams->min_segment_length = 0.0;

	// Add it to the culture
	AddNewSegment(seg);

	//restore old min seg length
	m_cultParams->min_segment_length = prev_min_seg_length;
}

//-----------------------------------------------------------------------------
// Create a new sprout force component for a newly formed branch
// TODO: What if both tips are active?
void Culture::CreateBranchingForce(Segment& seg)
{
	m_angio.total_bdyf = 0;																// Obtain the total number of sprouts
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
	recents.clear();
	//adding zero length segments shoudl be avoided
	//this will clear and refill recents
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
	m_nsegs++;
	m_frag.push_front(seg);
	assert(m_nsegs == (int)m_frag.size());
	//add the segment that was just created
	recents.push_back(&m_frag.front());

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
//FragmentBranching and Children from here

std::vector<FragmentBranching *> FragmentBranching::fragment_branchers;
std::set<FragmentBranching::BranchPoint> FragmentBranching::branch_points;
std::set<FragmentBranching::BranchPoint, FragmentBranching::BranchPointEpochCompare> FragmentBranching::parentgen;
bool FragmentBranching::create_timeline = true;
std::multiset<FragmentBranching::BranchPoint, FragmentBranching::BranchPointTimeFloorCompare> FragmentBranching::timeline;

//this is the combined growth step for all cultures/fragments/FragmentBranchers
void FragmentBranching::Grow()
{
	//get all segments which need values recalculated for them in this step
	if (!fragment_branchers.size())
		return;
	//Grow the segments populate the sets 
	for (size_t i = 0; i < fragment_branchers.size(); i++)
	{
		assert(fragment_branchers[i] != nullptr);
		fragment_branchers[i]->GrowSegments();
	}
	//keep track of the start iterators so that I can delete the processed data at the end of the step
	auto start_parentgen = parentgen.begin();

	auto start_branch_points = branch_points.begin();
	//generate all of the needed random numbers
	//note that branches length to branch will be generated at their emerge time
	SimulationTime end_time = fragment_branchers[0]->culture->m_pmat->m_pangio->CurrentSimTime();
	//the time in SimulationTime is the end of the current timestep
	assert((start_parentgen == parentgen.end()) || (start_parentgen->epoch_time <= end_time.t)); //parentgen should empty each timestep/ not be filled beyond the current step

	auto pg_iter = start_parentgen;
	auto br_iter = start_branch_points;

	while (true)
	{
		//kick out once all the points in this timestep have been processed
		if (pg_iter == parentgen.end() && br_iter == branch_points.end())
			break;
		if (pg_iter == parentgen.end() && br_iter->emerge_time > end_time.t)
			break;
		//find which iterator happened first if there  is a tie go one direction all the time
		if (br_iter == branch_points.end() || ((br_iter != branch_points.end() && pg_iter != parentgen.end()) &&((br_iter->emerge_time >= pg_iter->epoch_time) || ( br_iter->emerge_time > end_time.t))))
		{
			//go with pg_iter
			if (create_timeline)
			{
				BranchPoint pt = *pg_iter;
				pt.branch = false;
				timeline.insert(pt);
			}
			pg_iter->brancher->UpdateSegmentBranchDistance(pg_iter);
			++pg_iter;
			
		}
		else
		{
			//go with br iterator

			if (create_timeline)
			{
				BranchPoint pt = *br_iter;
				pt.branch = true;
				timeline.insert(pt);
			}
			//GrowSegment needs to populate the segment with branch times and recycle it if needed
			//copying the iterator protects this internal iterator
			br_iter->brancher->GrowSegment(br_iter);
			++br_iter;
		}
	}

	//delete the processed portion of the static structures should just reduce memory usage
	branch_points.erase(start_branch_points, br_iter);
	parentgen.erase(start_parentgen, pg_iter);
	//update the active tip lists
	for (size_t i = 0; i < fragment_branchers.size(); i++)
	{
		fragment_branchers[i]->culture->FindActiveTips();
	}
}



//does the portion of the grow that does not require rng
void FragmentBranching::GrowSegments()
{
	//update the vessel length and then grow the segments
	SimulationTime & st = culture->m_pmat->m_pangio->CurrentSimTime();
	const SegmentTipList & tiplist = culture->GetActiveTipList();
	auto iter = tiplist.begin();
	while (iter != tiplist.end())
	{
		GrowSegment(*iter, st.t - st.dt, st.dt);
		++iter;
	}

}

void NoFragmentBranching::GrowSegment(Segment::TIP * tip, double starttime, double grow_time)
{
	double dw[2];
	double da_value = culture->m_pmat->m_pangio->genericProjectToPoint(&tip->pt.ndomain->ElementRef(tip->pt.elemindex), &FEAngioNodeData::m_da, tip->pt.q);
	DirectionalWeights(da_value, dw);
	//TODO: this overwrite may not be ideal
	culture->m_pmat->m_cultureParams.vessel_orient_weights.x = dw[0];
	culture->m_pmat->m_cultureParams.vessel_orient_weights.y = dw[1];

	Segment seg = culture->GrowSegment(*tip, starttime, grow_time);
	culture->AddNewSegment(seg);
}

void ForwardFragmentBranching::GrowSegment(Segment::TIP * tip, double starttime, double grow_time)
{
	static int nseg_add = 0;
	double dw[2];
	double da_value = culture->m_pmat->m_pangio->genericProjectToPoint(&tip->pt.ndomain->ElementRef(tip->pt.elemindex), &FEAngioNodeData::m_da, tip->pt.q);
	DirectionalWeights(da_value, dw);
	//TODO: this overwrite may not be ideal
	culture->m_pmat->m_cultureParams.vessel_orient_weights.x = dw[0];
	culture->m_pmat->m_cultureParams.vessel_orient_weights.y = dw[1];

	Segment seg = culture->GrowSegment(*tip, starttime, grow_time);
	//now calculate the length to branch from the new tip
	seg.tip(1).length_to_branch = tip->length_to_branch;
	//unfortunatly the boundary condition code needs to add the branch/parent points as it may reposition branch points
	culture->AddNewSegment(seg);
	auto rseg = culture->RecentSegments();
	if (rseg.size() == 1) //handle segments in non bouncing mode
	{
		rseg[0]->tip(1).length_to_branch -= seg.length();
		if (rseg[0]->tip(1).length_to_branch < 0.0)
		{
			double bf = (rseg[0]->length() + rseg[0]->tip(1).length_to_branch) / rseg[0]->length();
			assert(bf >= 0.0 && bf <= 1.0);
			SimulationTime end_time = culture->m_pmat->m_pangio->CurrentSimTime();
			//TODO: this is too generous with the timing of the segment, consider adding time of death to segments
			double bt = mix(end_time.t - end_time.dt, end_time.t, bf);
			//add the points
			BranchPoint bp(bt, bt, rseg[0], bf, rseg[0]->m_nid, this);
			branch_points.insert(bp);
			parentgen.insert(bp);
		}
	}
	else if (rseg.size() == 0)
	{
		nseg_add++;
	}
}

void ForwardFragmentBranching::PostProcess(Segment & seg)
{
	//on initialization do nothing otherwise the segment is growing from 0 -> 1
	//update the time of birth of the segment
	
	if (!(seg.tip(0).bactive && seg.tip(1).bactive))
	{
		seg.SetTimeOfBirth(current_time);
		seg.tip(1).length_to_branch -= seg.length();
		if (seg.tip(1).length_to_branch < 0.0)
		{
			double bf = (seg.length() + seg.tip(1).length_to_branch) / seg.length();
			assert(bf >= 0.0 && bf <= 1.0);
			SimulationTime end_time = culture->m_pmat->m_pangio->CurrentSimTime();
			//TODO: this is too generous with the timing of the segment, consider adding time of death to segments
			double bt = mix(seg.GetTimeOfBirth(), end_time.t, bf);
			//add the points
			BranchPoint bp(bt, bt, &seg, bf, seg.m_nid, this);
			branch_points.insert(bp);
			parentgen.insert(bp);
		}
	}
	else
	{
		seg.SetTimeOfBirth(-1.0);//set the time of birth for initial segmnets to -1 
	}
	
}

void ForwardFragmentBranching::GrowSegment(std::set<BranchPoint>::iterator bp)
{
	//create a tip of origin
	static int mis_count = 0;
	Segment::TIP tip0 = bp->parent->tip(0);
	Segment::TIP tip1 = bp->parent->tip(1);
	vec3d pos = mix(bp->parent->tip(0).pos(), bp->parent->tip(1).pos(), bp->percent_of_parent);
	SimulationTime & st = culture->m_pmat->m_pangio->CurrentSimTime();
	if (culture->m_pmat->FindGridPoint(pos, tip0.pt.ndomain, tip0.pt.elemindex, tip0.pt))
	{
		//roll for the length to branch of the new segment
		//the overwrite here is okay as it is a copy
		tip0.length_to_branch = length_to_branch_point(culture->m_pmat->m_pangio->rengine);
		culture->BranchSegment(tip0, bp->emerge_time, st.t - bp->emerge_time);
	}
	else if (culture->m_pmat->FindGridPoint(pos, tip1.pt.ndomain, tip1.pt.elemindex, tip1.pt))
	{
		tip1.length_to_branch = length_to_branch_point(culture->m_pmat->m_pangio->rengine);
		culture->BranchSegment(tip1, bp->emerge_time, st.t - bp->emerge_time);
	}
	else
	{
		mis_count++;
		//consider failing when this is not in one of the elements containing the endpoints of the segments
		//assert(false);
	}
	if (culture->RecentSegments().size() == 0)
	{
		mis_count++;
	}
	
}
void ForwardFragmentBranching::UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp)
{
	double old_l2b = bp->parent->tip(1).length_to_branch;
	bp->parent->tip(1).length_to_branch += length_to_branch_point(culture->m_pmat->m_pangio->rengine);
	if (bp->parent->tip(1).length_to_branch < 0.0)
	{
		//add it back to the queue
		SimulationTime end_time = culture->m_pmat->m_pangio->CurrentSimTime();
		double bt = mix(bp->epoch_time, end_time.t, (old_l2b - bp->parent->tip(1).length_to_branch) / old_l2b);
		double bpct = (bt - (end_time.t - end_time.dt)) / (-(end_time.t - end_time.dt) + end_time.t);
		BranchPoint bpt(bt, bt, bp->parent, bpct , bp->priority, this);
		parentgen.insert(bpt);
		//make sure adding both back were wrong
		branch_points.insert(bpt);
	}
}
double ForwardFragmentBranching::GetLengthToBranch()
{
	return length_to_branch_point(culture->m_pmat->m_pangio->rengine);
}


//begin implementations of grow direction modifiers
vec3d GradientGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch)
{
	//calculate the density gradinet if above the threshold set the grow direction
	std::vector<double> densities;
	FESolidElement * se = dynamic_cast<FESolidElement*>(&tip.pt.ndomain->ElementRef(tip.pt.elemindex));
	densities = culture->m_pmat->m_pangio->createVectorOfMaterialParameters(se, &FEAngioNodeData::m_ecm_den);
	vec3d gradient = culture->m_pmat->m_pangio->gradient(se, densities, tip.pt.q);
	double gradnorm = gradient.norm();
	Segment seg;
	if (gradnorm > culture->m_pmat->m_cultureParams.density_gradient_threshold)
	{
		vec3d currentDirection = previous_dir;
		currentDirection.unit();
		vec3d currentDirectionGradientPlane = gradient ^ currentDirection;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ gradient;
		perpendicularToGradient.unit();
		return perpendicularToGradient;
	}
	return previous_dir;
}

vec3d AnastamosisGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch)
{
	return previous_dir;
}

vec3d BranchGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch)
{
	// If new segment is a branch we modify the grow direction a bit
	if (branch)
	{
		// TODO: what's the logic here? Why the 0.5 factor?
		//      If the vessel is aligned with the collagen (and the initial fragments are)
		//      then  the new branch will overlap the old segment.
		vec3d seg_vec = -previous_dir;
		vec3d coll_fib = culture->m_pmat->CollagenDirection(tip.pt);
		seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
		seg_vec.unit();
		return seg_vec;
	}
	return previous_dir;
}

vec3d DefaultGrowDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch)
{
	// Find the component of the new vessel direction determined by collagen fiber orientation    
	vec3d coll_dir = culture->m_pmat->CollagenDirection(tip.pt);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = mix(per_dir, coll_dir, culture->m_pmat->m_cultureParams.GetWeightInterpolation(culture->m_pmat->m_pangio->CurrentSimTime().dt));
	new_dir.unit();

	return new_dir;
}
