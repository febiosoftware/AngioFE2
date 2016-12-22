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
		/* removed until fully implemeented
	case 1:
		fbrancher = new ForwardFragmentBranching(this);
		break;
	case 2:
		fbrancher = new PsuedoDeferedFragmentBranching(this);
		break;
		*/
	default:
		assert(false);
	}

	// do the initial seeding
	if (!fseeder->SeedFragments(m_angio.CurrentSimTime(), this))
		return false;
	previous_average_segment_length = m_cultParams->average_length_to_branch_point;
	previous_standard_deviation = m_cultParams->std_deviation;
	segment_length_distribution = normal_distribution<double>(previous_average_segment_length, previous_standard_deviation);

	FragmentBranching::SetupBranchPoints();
	FragmentBranching::GrowInitialBranches();
	return true;
}
void Culture::SetLengthToBranch()
{
	auto iter = m_frag.begin();
	while (iter != m_frag.end())
	{
		//note the lenghts are divided by two so the average branch point distance is maintained
		//TODO: this maintains the average length to branch for the entire vessel but the standard deviation of initial branch points is halved
		//shoudl this use a custom distribution for seeding
		iter->tip(0).length_to_branch = (segment_length_distribution(m_pmat->m_pangio->rengine) - iter->length()) / 2;
		iter->tip(1).length_to_branch = (segment_length_distribution(m_pmat->m_pangio->rengine) - iter->length()) / 2;
		++iter;
	}
}

void Culture::InitialBranching()
{
	//consider that the tips are growing, this creates a special case when seeding the fragments
	//the tip does matter as both tips could be up for regrowth
	std::map<double, std::map<int, Segment::TIP*>> regrowth;
	bool new_branch = false;
	//segments must be processed in chronological order with ties being broken by vessel_id
	auto iter = m_frag.begin();

	while (iter != m_frag.end())
	{
		for (int i = 0; i < 2; i++)
		{
			//if the tip is not active do nothing
			if ((iter->tip(i).length_to_branch < 0.0) && iter->tip(i).bactive)
			{
				//consider making the following a function
				//calculate the ratio of time when this segment diveged
				double bt = ((iter->length()/2) + iter->tip(i).length_to_branch) / (iter->length()/2);

				//values closer to zero happened earlier chronologically
				//make sure that the branch should have happened in this grow step
				assert(bt >= 0 && bt <= 1.0);
				if (regrowth.count(bt))
				{
					regrowth.find(bt)->second.insert(std::pair<int,  Segment::TIP*>(iter->m_nid, &iter->tip(i)));
				}
				else
				{
					std::map<int, Segment::TIP*> temp;
					regrowth.insert(std::pair<double, std::map<int, Segment::TIP*>>(bt, temp));
					regrowth.find(bt)->second.insert(std::pair<int, Segment::TIP*>(iter->m_nid, &(iter->tip(i))));
				}
			}
		}
		++iter;
	}

	//now create branches and recalculate distance to branch
	auto iterr = regrowth.begin();
	while (iterr != regrowth.end())
	{
		auto iterk = iterr->second.begin();
		while (iterk != iterr->second.end())
		{
			//branch length
			double bl = (1.0 - iterr->first)* iterk->second->parent->length()/2;
			assert(bl > 0.0);
			//add the branch 
			Segment nseg = *iterk->second->parent;
			//update the branch point of the base segment
			iterk->second->length_to_branch += segment_length_distribution(m_pmat->m_pangio->rengine);
			if (iterk->second->length_to_branch < 0.0)
			{
				//add this back to regrowth
				//calculate the ratio of time when this segment diveged
				double bt = ((iterk->second->parent->length()/2) + iterk->second->length_to_branch) / (iterk->second->parent->length()/2);

				//values closer to zero happened earlier chronologically
				//make sure that the branch should have happened in this grow step
				assert(bt >= 0 && bt <= 1.0);
				if (regrowth.count(bt))
				{
					regrowth.find(bt)->second.insert(std::pair<int, Segment::TIP*>(iterk->second->parent->m_nid, iterk->second));
				}
				else
				{
					std::map<int, Segment::TIP*> temp;
					regrowth.insert(std::pair<double, std::map<int, Segment::TIP*>>(bt, temp));
					regrowth.find(bt)->second.insert(std::pair<int, Segment::TIP*>(iterk->second->parent->m_nid, iterk->second));
				}
			}

			//find the position of the astil
			//check the order of the tips
			//consider only the portion of the segment from the midpoint to the current tip
			vec3d astil = mix(mix(iterk->second->parent->tip(0).pt.r, iterk->second->parent->tip(1).pt.r, 0.5), iterk->second->pos(), iterr->first);
			Segment::TIP nt0 = iterk->second->parent->tip(0);
			Segment::TIP nt1 = iterk->second->parent->tip(1);

			
			
			if (m_pmat->FindGridPoint(astil, nt0.pt.ndomain, nt0.pt.elemindex, nt0.pt))
			{
				nseg = GrowSegmentBranching(nt0, bl);
				AddNewSegment(nseg);
			}
			else if (m_pmat->FindGridPoint(astil, nt1.pt.ndomain, nt1.pt.elemindex, nt1.pt))
			{
				nseg = GrowSegmentBranching(nt1, bl);
				AddNewSegment(nseg);
			}
			else
			{
				assert(false);
			}
			//increment the branch count
			m_num_branches++;

			//calculate the distance until the branch branches
			double utnb = segment_length_distribution(m_pmat->m_pangio->rengine);

			//if the segment would go out of bounds it needs special handling
			if (nseg.GetFlag(Segment::BC_DEAD))
			{
				//special handling of encountering a boundary condition
				//nothing needs done on stop
				//bouncy and passthrough will be more difficult
			}
			else
			{
				nseg.tip(1).length_to_branch += utnb;
				if (nseg.tip(1).length_to_branch < 0.0)
				{
					//add this back to the datastructure

					double bt = (nseg.length() + nseg.tip(1).length_to_branch) / nseg.length();

					//values closer to zero happened earlier chronologically
					//make sure that the branch should have happened in this grow step
					assert(bt >= 0 && bt <= 1.0);
					if (regrowth.count(bt))
					{
						regrowth.find(bt)->second.insert(std::pair<int, Segment::TIP *>(iter->m_nid, &nseg.tip(1)));
					}
					else
					{
						std::map<int, Segment::TIP *> temp;
						regrowth.insert(std::pair<double, std::map<int, Segment::TIP *>>(bt, temp));
						regrowth.find(bt)->second.insert(std::pair<int, Segment::TIP *>(iter->m_nid, &nseg.tip(1)));
					}
				}
			}


			++iterk;
		}
		
		++iterr;
	}

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
		vec3d q = vrand();

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

	// decide if this initial segment is allowed to branch
	if (frand() < culture_params->m_initial_branch_probability) seg.SetFlagOn(Segment::INIT_BRANCH);

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
//updates the segment length distribution if needed
void Culture::UpdateSegmentDistribution()
{
	if ((previous_average_segment_length != m_cultParams->average_length_to_branch_point) || (previous_standard_deviation != m_cultParams->std_deviation))
	{
		//set the new distributio and previous values
		//consider setting a breakpoint here as this changing at the wrong time could give diverging networks
		//TODO: this may need revised if another distribution is used
		segment_length_distribution = normal_distribution<double>(m_cultParams->average_length_to_branch_point, m_cultParams->std_deviation);
		previous_average_segment_length = m_cultParams->average_length_to_branch_point;
		previous_standard_deviation = m_cultParams->std_deviation;
	}
}

//-----------------------------------------------------------------------------
// Grow phase
void Culture::GrowVessels()
{
	FEMesh* mesh = m_angio.GetMesh();

	for (TipIter it = m_active_tips.begin(); it != m_active_tips.end(); ++it)
	{
		Segment::TIP& tip = *(*it);
		if (!tip.bactive)
			continue; //this may happen when the segment has grown into the boundary but has not yet been removed from the active tip list
		//modify direction based on anisotropy value
		double dw[2];
		double da_value = m_angio.genericProjectToPoint(&tip.pt.ndomain->ElementRef(tip.pt.elemindex), &FEAngioNodeData::m_da, tip.pt.q);
		DirectionalWeights(da_value, dw);
		//TODO: this overwrite may not be ideal
		m_cultParams->vessel_orient_weights.x = dw[0];
		m_cultParams->vessel_orient_weights.y = dw[1];
		
		//calculate the density gradinet if above the threshold set the grow direction
		std::vector<double> densities;
		FESolidElement * se = dynamic_cast<FESolidElement*>(&tip.pt.ndomain->ElementRef(tip.pt.elemindex));
		densities = m_angio.createVectorOfMaterialParameters(se, &FEAngioNodeData::m_ecm_den);
		vec3d gradient = m_angio.gradient(se, densities, tip.pt.q);
		double gradnorm = gradient.norm();
		Segment seg;
		if (gradnorm > m_cultParams->density_gradient_threshold)
		{
			vec3d currentDirection = FindGrowDirection(tip);
			currentDirection.unit();
			vec3d currentDirectionGradientPlane = gradient ^ currentDirection;
			currentDirectionGradientPlane.unit();
			vec3d perpendicularToGradient = currentDirectionGradientPlane ^ gradient;
			perpendicularToGradient.unit();
			seg = GrowSegment(tip, false, false, perpendicularToGradient);
			AddNewSegment(seg);
		}
		else
		{
			// Create new vessel segment at the current tip existing segment
			seg = GrowSegment(tip);
			AddNewSegment(seg);
		}
		

		// Add the segment to the network
		// This will also enforce the boundary conditions
		
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
		vec3d coll_fib = m_pmat->CollagenDirection(tip.pt);
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

Segment Culture::GrowSegmentBranching(Segment::TIP& tip, double length)
{
	// this is the new segment length
	double seg_length = length;

	// determine the growth direction
	vec3d seg_vec = FindGrowDirection(tip);

	//      If the vessel is aligned with the collagen (and the initial fragments are)
	//      then  the new branch will overlap the old segment. 
	vec3d coll_fib = m_pmat->CollagenDirection(tip.pt);
	seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
	seg_vec.unit();


	// Create a new segment
	Segment seg;

	// transer seed label
	seg.seed(tip.nseed);

	seg.vessel(m_num_vessel++);

	// create a new vessel
	// the first point is a copy of the last point of the source segment
	seg.tip(0).pt = tip.pt;

	// position the end point
	seg.tip(1).pt.r = seg.tip(0).pos() + seg_vec*seg_length;

	seg.tip(0).bactive = false;		// Turn off origin tip of new segment
	seg.tip(1).bactive = true;		// Turn on end tip of new segment

	// Turn off previous segment tip
	tip.bactive = false;

	m_pmat->FindGridPoint(seg.tip(1).pos(), seg.tip(1).pt);

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
			double bprob = den_scale*m_cultParams->GetBranchProbility(time.dt);
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
// Branching phase of the growth step.
void Culture::BranchVessels2(SimulationTime& time)
{
	std::map<double, std::map<int, Segment *> > remaining_segments;//uses two keys the time and the vessel id
	// Elongate the active vessel tips
	for (SegIter it = m_frag.begin(); it != m_frag.end(); ++it)
	{
		// decide if this segment can branch
		bool branch = true;
		if (it->GetFlag(Segment::BC_DEAD)) branch = false;	// Make sure segment has not reached a boundary
		if (it->GetFlag(Segment::ANAST)) branch = false;	// and has not undegone anastimoses.
		if (it->tip(0).bactive || it->tip(1).bactive) branch = false;	// don't branch segments that are actively growing
		int vessel = it->m_nvessel;
		if (branch || it->GetFlag(Segment::INIT_BRANCH))
		{
			// pick an end randomly
			int k = (frand() < 0.5 ? k = 0 : k = 1);

			// find the density scale factor
			double den_scale = FindDensityScale(it->tip(k).pt);

			// calculate actual branching probability
			double bprob = den_scale*m_cultParams->GetBranchProbility(time.dt);
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

	//calculate the density gradinet if above the threshold set the grow direction
	std::vector<double> densities;
	FESolidElement * se = dynamic_cast<FESolidElement*>(&tip.pt.ndomain->ElementRef(tip.pt.elemindex));
	densities = m_angio.createVectorOfMaterialParameters(se, &FEAngioNodeData::m_ecm_den);
	vec3d gradient = m_angio.gradient(se, densities, tip.pt.q);
	double gradnorm = gradient.norm();
	Segment seg;
	if (gradnorm > m_cultParams->density_gradient_threshold)
	{
		vec3d currentDirection = FindGrowDirection(tip);
		currentDirection.unit();
		vec3d currentDirectionGradientPlane = gradient ^ currentDirection;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ gradient;
		perpendicularToGradient.unit();
		//do a special calculation 
		vec3d coll_fib = m_pmat->CollagenDirection(tip.pt);
		vec3d seg_vec = perpendicularToGradient;
		seg_vec = coll_fib - seg_vec*(seg_vec*coll_fib)*0.5;
		seg_vec.unit();
		gradient.unit();
		//project the direction onto the plane u - proj_v(u)
		//note the mag of gradeint is 1
		seg_vec = seg_vec - (gradient * (seg_vec * gradient));

		seg = GrowSegment(tip, false, true, seg_vec);
	}
	else
	{
		// Create new vessel segment at the current tip existing segment
		seg = GrowSegment(tip, true, true);
	}

	// Add it to the culture
	AddNewSegment(seg);
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
    vec3d coll_dir = m_pmat->CollagenDirection(tip.pt);

	// Component of new vessel orientation resulting from previous vessel direction        
	vec3d per_dir = tip.u;

	vec3d new_dir = mix(per_dir, coll_dir, m_cultParams->GetWeightInterpolation(m_pmat->m_pangio->CurrentSimTime().dt));
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
//FragmentBranching and Children from here

std::vector<FragmentBranching *> FragmentBranching::fragment_branchers;
std::set<FragmentBranching::BranchPoint> FragmentBranching::branch_points;
std::set<FragmentBranching::BranchPoint, FragmentBranching::BranchPointEpochCompare> FragmentBranching::parentgen;

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
	assert((start_parentgen == parentgen.end()) || (start_parentgen->epoch_time < end_time.t)); //parentgen should empty each timestep/ not be filled beyond the current step

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
		if (pg_iter == parentgen.end() || ((br_iter != branch_points.end()) &&(br_iter->emerge_time < pg_iter->epoch_time)))
		{
			//go with br iterator

			//GrowSegment needs to populate the segment with branch times and recycle it if needed
			//copying the iterator protects this internal iterator
			br_iter->brancher->GrowSegment(br_iter);
			++br_iter;
		}
		else if (br_iter == branch_points.end() || br_iter->emerge_time > end_time.t)
		{
			//go with pg_iter
			pg_iter->brancher->UpdateSegmentBranchDistance(pg_iter);
			++pg_iter;
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

void FragmentBranching::GrowInitialBranches()
{
	//same as grow except for there is no initial growth of all active tips just the cleanup of random numbers

	//get all segments which need values recalculated for them in this step
	if (!fragment_branchers.size())
		return;
	//keep track of the start iterators so that I can delete the processed data at the end of the step
	auto start_parentgen = parentgen.begin();

	auto start_branch_points = branch_points.begin();
	//generate all of the needed random numbers
	//note that branches length to branch will be generated at their emerge time
	SimulationTime end_time = fragment_branchers[0]->culture->m_pmat->m_pangio->CurrentSimTime();
	//the time in SimulationTime is the end of the current timestep
	assert((start_parentgen == parentgen.end()) || (start_parentgen->epoch_time < end_time.t)); //parentgen should empty each timestep/ not be filled beyond the current step

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
		if (pg_iter == parentgen.end() || ((br_iter != branch_points.end()) && (br_iter->emerge_time < pg_iter->epoch_time)))
		{
			//go with br iterator

			//GrowSegment needs to populate the segment with branch times and recycle it if needed
			//copying the iterator protects this internal iterator
			br_iter->brancher->GrowSegment(br_iter);
			++br_iter;
		}
		else if (br_iter == branch_points.end() || br_iter->emerge_time > end_time.t)
		{
			//go with pg_iter
			pg_iter->brancher->UpdateSegmentBranchDistance(pg_iter);
			++pg_iter;
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
	culture->UpdateNewVesselLength(culture->m_pmat->m_pangio->CurrentSimTime());
	const SegmentTipList & tiplist = culture->GetActiveTipList();
	auto iter = tiplist.begin();
	while (iter != tiplist.end())
	{
		GrowSegment(*iter);
		++iter;
	}

}

void NoFragmentBranching::GrowSegment(Segment::TIP * tip)
{
	double dw[2];
	double da_value = culture->m_pmat->m_pangio->genericProjectToPoint(&tip->pt.ndomain->ElementRef(tip->pt.elemindex), &FEAngioNodeData::m_da, tip->pt.q);
	DirectionalWeights(da_value, dw);
	//TODO: this overwrite may not be ideal
	culture->m_pmat->m_cultureParams.vessel_orient_weights.x = dw[0];
	culture->m_pmat->m_cultureParams.vessel_orient_weights.y = dw[1];

	//calculate the density gradinet if above the threshold set the grow direction
	std::vector<double> densities;
	FESolidElement * se = dynamic_cast<FESolidElement*>(&tip->pt.ndomain->ElementRef(tip->pt.elemindex));
	densities = culture->m_pmat->m_pangio->createVectorOfMaterialParameters(se, &FEAngioNodeData::m_ecm_den);
	vec3d gradient = culture->m_pmat->m_pangio->gradient(se, densities, tip->pt.q);
	double gradnorm = gradient.norm();
	Segment seg;
	if (gradnorm > culture->m_pmat->m_cultureParams.density_gradient_threshold)
	{
		vec3d currentDirection = culture->FindGrowDirection(*tip);
		currentDirection.unit();
		vec3d currentDirectionGradientPlane = gradient ^ currentDirection;
		currentDirectionGradientPlane.unit();
		vec3d perpendicularToGradient = currentDirectionGradientPlane ^ gradient;
		perpendicularToGradient.unit();
		seg = culture->GrowSegment(*tip, false, false, perpendicularToGradient);
		culture->AddNewSegment(seg);
	}
	else
	{
		// Create new vessel segment at the current tip existing segment
		seg = culture->GrowSegment(*tip);
		culture->AddNewSegment(seg);
	}
}

void FragmentBranching::SetupBranchPoints()
{
	//Grow the segments populate the sets 
	for (size_t i = 0; i < fragment_branchers.size(); i++)
	{
		auto tiplist = fragment_branchers[i]->culture->GetActiveTipList();
		auto iter = tiplist.begin();
		while (iter != tiplist.end())
		{
			fragment_branchers[i]->InitializeSegment(*iter);
			++iter;
		}
	}
	
}