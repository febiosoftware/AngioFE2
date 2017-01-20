#include "FragmentBranching.h"
#include "angio3d.h"
#include "FEAngio.h"

//FragmentBranching and Children from here

std::vector<FragmentBranching *> FragmentBranching::fragment_branchers;
std::set<FragmentBranching::BranchPoint> FragmentBranching::branch_points;
std::set<FragmentBranching::BranchPoint, FragmentBranching::BranchPointEpochCompare> FragmentBranching::parentgen;
bool FragmentBranching::create_timeline = true;
std::multiset<FragmentBranching::BranchPoint, FragmentBranching::BranchPointTimeFloorCompare> FragmentBranching::timeline;

void FragmentBranching::SetCulture(Culture * cp)
{
	culture = cp;
}

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
		if (br_iter == branch_points.end() || ((br_iter != branch_points.end() && pg_iter != parentgen.end()) && ((br_iter->emerge_time >= pg_iter->epoch_time) || (br_iter->emerge_time > end_time.t))))
		{
			//go with pg_iter
#ifndef NDEBUG
				BranchPoint pt = *pg_iter;
				pt.branch = false;
				timeline.insert(pt);
#endif
			pg_iter->brancher->UpdateSegmentBranchDistance(pg_iter);
			++pg_iter;

		}
		else
		{
			//go with br iterator

#ifndef NDEBUG
			BranchPoint pt = *br_iter;
			pt.branch = true;
			timeline.insert(pt);
#endif
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
//calculates the time a segment took to grow during the current timestep, this assumes the growthrate is constant within a step
double FragmentBranching::TimeOfGrowth(Segment * seg)
{
	
	double avg_den = culture->m_pmat->m_pangio->FindECMDensity(seg->tip(0).pt);//mix(culture->m_pmat->m_pangio->FindECMDensity(seg->tip(0).pt), culture->m_pmat->m_pangio->FindECMDensity(seg->tip(1).pt), 0.5);
	double den_scale = culture->m_pmat->m_cultureParams.m_density_scale_factor.x + culture->m_pmat->m_cultureParams.m_density_scale_factor.y
		*exp(-culture->m_pmat->m_cultureParams.m_density_scale_factor.z*avg_den);
	SimulationTime & time = culture->m_pmat->m_pangio->CurrentSimTime();
	double base_length = culture->SegmentLength(time.t - time.dt, time.dt);
	double rv = time.dt *( seg->length() / (den_scale *base_length));
	assert(rv > 0.0);
	assert(rv <= (1.01*time.dt));
	return rv;
}

void DirectionalWeights(double da, double dw[2])
{
	double degree_isotropy = 1.0 - da;
	dw[0] = degree_isotropy;
	dw[1] = da;
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




//begin implementation of psuedo deferred branching
PsuedoDeferedFragmentBranching::PsuedoDeferedFragmentBranching(FEModel * model) : FragmentBranching(model)
{
	AddProperty(&length_to_branch_point, "length_to_branch");
	AddProperty(&time_to_emerge, "time_to_emerge");
}
void PsuedoDeferedFragmentBranching::GrowSegment(Segment::TIP * tip, double starttime, double grow_time)
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
	seg.tip(1).wait_time_to_branch = tip->wait_time_to_branch;
	culture->AddNewSegment(seg);
	auto rseg = culture->RecentSegments();
	if (rseg.size() == 1) //handle segments in non bouncing mode
	{
		rseg[0]->tip(1).length_to_branch -= seg.length();
		rseg[0]->SetTimeOfBirth(starttime);
		if (rseg[0]->tip(1).length_to_branch < 0.0)
		{
			double bf = (rseg[0]->length() + rseg[0]->tip(1).length_to_branch) / rseg[0]->length();
			assert(bf >= 0.0 && bf <= 1.0);
			SimulationTime end_time = culture->m_pmat->m_pangio->CurrentSimTime();
			double bt = mix(end_time.t - end_time.dt, end_time.t - end_time.dt + TimeOfGrowth(rseg[0]), bf);
			//add the points
			assert(tip->wait_time_to_branch >= 0.0);
			BranchPoint bp(bt + tip->wait_time_to_branch, bt, rseg[0], bf, rseg[0]->m_nid, this, 4);
			branch_points.insert(bp);
			parentgen.insert(bp);
		}
	}
	else if (rseg.size() == 0)
	{
		nseg_add++;
	}
	else
	{

	}
}

void PsuedoDeferedFragmentBranching::PostProcess(Segment & seg)
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
			BranchPoint bp(bt + seg.tip(1).wait_time_to_branch, bt, &seg, bf, seg.m_nid, this, 5);
			branch_points.insert(bp);
			parentgen.insert(bp);
		}
	}
	else
	{
		seg.SetTimeOfBirth(-1.0);//set the time of birth for initial segmnets to -1 
	}

}

void PsuedoDeferedFragmentBranching::GrowSegment(std::set<BranchPoint>::iterator bp)
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
		tip0.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
		tip0.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
		culture->BranchSegment(tip0, bp->emerge_time, st.t - bp->emerge_time);
	}
	else if (culture->m_pmat->FindGridPoint(pos, tip1.pt.ndomain, tip1.pt.elemindex, tip1.pt))
	{
		tip1.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
		tip1.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
		culture->BranchSegment(tip1, bp->emerge_time, st.t - bp->emerge_time);
	}
	else
	{
		//this is when the segment is in another element ie the segment is growing through at least 2 elements
		//this may indicate a problem in getting consistent results
		if (culture->m_pmat->FindGridPoint(pos, tip0.pt))
		{
			tip0.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
			tip0.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
			culture->BranchSegment(tip0, bp->emerge_time, st.t - bp->emerge_time);

		}
		else
		{
			//failed to place the tip in another element
			assert(false);
		}
	}
	ProcessNewSegments(bp->emerge_time);

}

void PsuedoDeferedFragmentBranching::ProcessNewSegments(double start_time)
{
	auto rseg = culture->RecentSegments();

	//consider handling when no segments are added
	if (rseg.size() == 0)
	{

	}
	//the easy case just a single segment has been added
	else
	{
		double ct = start_time;
		for (size_t i = 0; i < rseg.size(); i++)
		{
			rseg[i]->SetTimeOfBirth(ct);
			ct += TimeOfGrowth(rseg[i]);
		}
	}
}
void PsuedoDeferedFragmentBranching::UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp)
{
	double old_l2b = bp->parent->tip(1).length_to_branch;
	bp->parent->tip(1).length_to_branch += length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
	if (bp->parent->tip(1).length_to_branch < 0.0)
	{
		//add it back to the queue
		SimulationTime end_time = culture->m_pmat->m_pangio->CurrentSimTime();
		double bt = mix(bp->epoch_time, end_time.t, (old_l2b - bp->parent->tip(1).length_to_branch) / old_l2b);
		double bpct = (bt - (end_time.t - end_time.dt)) / (-(end_time.t - end_time.dt) + end_time.t);
		BranchPoint bpt(bt + bp->parent->tip(1).wait_time_to_branch, bt, bp->parent, bpct, bp->priority, this, 6);
		parentgen.insert(bpt);
		//make sure adding both back were wrong
		branch_points.insert(bpt);
	}
}
double PsuedoDeferedFragmentBranching::GetLengthToBranch()
{
	double rv = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
	assert(rv >= 0.0);
	return rv;
}
double PsuedoDeferedFragmentBranching::GetTimeToEmerge()
{
	double rv = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
	assert(rv >= 0.0);
	return rv;
}