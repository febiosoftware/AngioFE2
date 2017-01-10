#include "FragmentBranching.h"
#include "angio3d.h"
#include "FEAngio.h"

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
		if (br_iter == branch_points.end() || ((br_iter != branch_points.end() && pg_iter != parentgen.end()) && ((br_iter->emerge_time >= pg_iter->epoch_time) || (br_iter->emerge_time > end_time.t))))
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
		BranchPoint bpt(bt, bt, bp->parent, bpct, bp->priority, this);
		parentgen.insert(bpt);
		//make sure adding both back were wrong
		branch_points.insert(bpt);
	}
}
double ForwardFragmentBranching::GetLengthToBranch()
{
	return length_to_branch_point(culture->m_pmat->m_pangio->rengine);
}