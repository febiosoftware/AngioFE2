#include "StdAfx.h"
#include "FragmentSeeder.h"
#include "Segment.h"
#include "Culture.h"
#include "FEAngio.h"
#include "angio3d.h"
#include "FECore/FEMesh.h"



FragmentSeeder::FragmentSeeder(FEModel * model) : FEMaterial(model), culture(nullptr)
{

}
void FragmentSeeder::SetCulture(Culture * cp)
{
	culture = cp;
}

BEGIN_PARAMETER_LIST(FragmentSeeder, FEMaterial)
ADD_PARAMETER(number_fragments, FE_PARAM_INT, "number_fragments");
ADD_PARAMETER(initial_vessel_length, FE_PARAM_DOUBLE, "initial_vessel_length");
END_PARAMETER_LIST();
//-----------------------------------------------------------------------------
// Create initial fragments

ClassicFragmentSeeder::ClassicFragmentSeeder(FEModel * model) : FragmentSeeder(model)
{
	
}

bool ClassicFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	for (int i = 0; i < number_fragments; ++i)
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
	culture->m_num_vessel = number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	return true;
}

//-----------------------------------------------------------------------------
// Generates an initial fragment that lies inside the grid.
// the classic mode still uses a biased growth direction
bool ClassicFragmentSeeder::createInitFrag(Segment& seg)
{
	FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
	// Set seg length to value of growth function at t = 0
	double seg_length = initial_vessel_length;

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
		p0 = culture->m_pmat->m_pangio->FindGridPoint(dynamic_cast<FESolidDomain*>(&mesh->Domain(0)), elem_num, q);

		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = vrand();
		seg_vec.unit();

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		culture->m_pmat->FindGridPoint(r1, p1);

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

	// Mark segment as an initial fragment
	seg.SetFlagOn(Segment::INIT_SPROUT);

	// all good
	return true;
}


MultiDomainFragmentSeeder::MultiDomainFragmentSeeder(FEModel * model) : FragmentSeeder(model)
{

}
//-----------------------------------------------------------------------------
// Create initial fragments
bool MultiDomainFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
	int elementsInDomains = 0;
	for (size_t i = 0; i < culture->m_pmat->domains.size(); i++)
	{
		elementsInDomains += culture->m_pmat->m_pangio->GetMesh()->Domain(culture->m_pmat->domains[i]).Elements();
		domains.emplace_back(&culture->m_pmat->m_pangio->GetMesh()->Domain(culture->m_pmat->domains[i]));
	}

	std::uniform_int_distribution<int> ddist(0, culture->m_pmat->domains.size() - 1);
	std::uniform_int_distribution<int> edist(0, elementsInDomains - 1);
	SegGenItem sgi;

	for (int i = 0; i < number_fragments; ++i)
	{
		// Create an initial segment
		Segment seg;
		sgi.domain = dynamic_cast<FESolidDomain*>(&mesh->Domain(culture->m_pmat->domains[ddist(culture->m_pmat->m_pangio->rengine)]));
		sgi.ielement = edist(culture->m_pmat->m_pangio->rengine);
		if (createInitFrag(seg, sgi, culture) == false) return false;

		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		// add it to the list
		culture->AddSegment(seg);
	}

	// init vessel counter
	culture->m_num_vessel = number_fragments;

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
	double seg_length = initial_vessel_length;

	// Since the creation of such a segment may fail (i.e. too close to boundary)
	// we loop until we find one.
	const int MAX_TRIES = 10;
	int ntries = 0;
	GridPoint p0, p1;
	do
	{
		// generate random natural coordinates
		vec3d q = culture->m_pmat->m_pangio->uniformInUnitCube();

		// set the position of the first tip
		p0 = culture->m_pmat->m_pangio->FindGridPoint(item.domain, item.ielement, q);

		// Determine vessel orientation based off of collagen fiber orientation
		vec3d seg_vec = culture->m_pmat->m_pangio->uniformRandomDirection();
		seg_vec.unit();

		// End of new segment is origin plus length component in each direction	
		vec3d r1 = p0.r + seg_vec*seg_length;

		// find the element where the second tip is
		culture->m_pmat->FindGridPoint(r1, p1);
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
	seg.tip(0).wait_time_to_branch = culture->fbrancher->GetTimeToEmerge();
	assert(seg.tip(0).wait_time_to_branch > 0.0);
	seg.tip(1).wait_time_to_branch = culture->fbrancher->GetTimeToEmerge();
	assert(seg.tip(1).wait_time_to_branch > 0.0);

	assert(seg.tip(0).wait_time_to_branch >= 0);
	assert(seg.tip(1).wait_time_to_branch >= 0);

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
		return findElement(val, lo, mid - 1, begin, end);
	}
	else if (val > end[mid])
	{
		return findElement(val, mid + 1, high, begin, end);
	}
	else
	{
		return mid;
	}
}

MDByVolumeFragmentSeeder::MDByVolumeFragmentSeeder(FEModel * model) :FragmentSeeder(model)
{
	
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool MDByVolumeFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
	int elementsInDomains = 0;
	for (size_t i = 0; i < culture->m_pmat->domains.size(); i++)
	{
		elementsInDomains += culture->m_pmat->m_pangio->GetMesh()->Domain(culture->m_pmat->domains[i]).Elements();
		domains.emplace_back(&culture->m_pmat->m_pangio->GetMesh()->Domain(culture->m_pmat->domains[i]));
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

	for (int i = 0; i < number_fragments; ++i)
	{
		//do a binary search to find the element that contains the volume
		double vol = voluchoice(culture->m_pmat->m_pangio->rengine);
		size_t ei = findElement(vol, 0, elementsInDomains - 1, totalWeightsBegin, totalWeightsEnd);
		// Create an initial segment
		Segment seg;
		FEElement * elem = elements[ei];
		sgi.ielement = elem->GetID() - 1 - culture->m_pmat->meshOffsets[elem->GetDomain()];
		sgi.domain = dynamic_cast<FESolidDomain*>(elem->GetDomain());

		if (createInitFrag(seg, sgi, culture) == false) return false;

		// Give the segment the appropriate label
		seg.seed(i);

		// Set the segment vessel as the segment label
		seg.vessel(seg.seed());

		// add it to the list
		culture->AddSegment(seg);
	}

	// init vessel counter
	culture->m_num_vessel = number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	delete[] totalWeightsBegin;
	delete[] totalWeightsEnd;
	delete[] elements;
	return true;
}

MDAngVessFileFragmentSeeder::MDAngVessFileFragmentSeeder(FEModel * model) :FragmentSeeder(model)
{
	vessel_file[0] = 0;
}

//-----------------------------------------------------------------------------
// Create initial fragments
bool MDAngVessFileFragmentSeeder::SeedFragments(SimulationTime& time, Culture * culture)
{
	FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
	infile.open(vessel_file);
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
		int festate = 0;
		double length = 0.0f, segtime = 0;
		line = regex_replace(line, whitespace, " ");
		if (9 != sscanf(line.c_str(), "%d %lf %lf %lf %lf %lf %lf %lf %lf", &festate, &segtime, &p1.x, &p1.y, &p1.z, &p2.x, &p2.y, &p2.z, &length))
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
	culture->m_num_vessel = number_fragments;

	// Update the active tip container
	culture->FindActiveTips();

	return true;
}
BEGIN_PARAMETER_LIST(MDAngVessFileFragmentSeeder, FragmentSeeder)
ADD_PARAMETER(vessel_file, FE_PARAM_STRING, "vessel_file");
END_PARAMETER_LIST();