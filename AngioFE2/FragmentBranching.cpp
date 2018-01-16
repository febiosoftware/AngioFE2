#include "FragmentBranching.h"
#include "angio3d.h"
#include "FEAngio.h"
#include "FEAngioMaterial.h"

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
void FragmentBranching::Grow(double start_time, double grow_time)
{
	//get all segments which need values recalculated for them in this step
	if (!fragment_branchers.size())
		return;

	//Grow the segments populate the sets 
	for (size_t i = 0; i < fragment_branchers.size(); i++)
	{
		assert(fragment_branchers[i] != nullptr);
		fragment_branchers[i]->GrowSegments(start_time, grow_time);
	}
	//keep track of the start iterators so that I can delete the processed data at the end of the step
	auto start_parentgen = parentgen.begin();

	auto start_branch_points = branch_points.begin();
	//generate all of the needed random numbers
	//note that branches length to branch will be generated at their emerge time
	double end_time = start_time + grow_time;
	//the time in SimulationTime is the end of the current timestep
	assert((start_parentgen == parentgen.end()) || (start_parentgen->epoch_time <= end_time)); //parentgen should empty each timestep/ not be filled beyond the current step

	auto pg_iter = start_parentgen;
	auto br_iter = start_branch_points;

	while (true)
	{
		//kick out once all the points in this timestep have been processed
		if (pg_iter == parentgen.end() && br_iter == branch_points.end())
			break;
		if (pg_iter == parentgen.end() && br_iter->emerge_time > end_time)
			break;
		//find which iterator happened first if there  is a tie go one direction all the time
		if (br_iter == branch_points.end() || ((br_iter != branch_points.end() && pg_iter != parentgen.end()) && ((br_iter->emerge_time >= pg_iter->epoch_time) || (br_iter->emerge_time > end_time))))
		{
			//go with pg_iter
#ifndef NDEBUG
				BranchPoint pt = *pg_iter;
				timeline.insert(pt);
#endif
			pg_iter->brancher->UpdateSegmentBranchDistance(pg_iter, start_time, grow_time);
			++pg_iter;

		}
		else
		{
			//go with br iterator

#ifndef NDEBUG
			BranchPoint pt = *br_iter;
			timeline.insert(pt);
#endif
			//GrowSegment needs to populate the segment with branch times and recycle it if needed
			//copying the iterator protects this internal iterator
			br_iter->brancher->GrowSegment(br_iter, start_time, grow_time);
			++br_iter;
		}
	}

	//delete the processed portion of the static structures should just reduce memory usage
	branch_points.erase(start_branch_points, br_iter);
	parentgen.erase(start_parentgen, pg_iter);
	//update the active tip lists
	for (int i = 0; i < fragment_branchers.size(); i++)
	{
		fragment_branchers[i]->culture->FindActiveTips();
	}
}



//does the portion of the grow that does not require rng
void FragmentBranching::GrowSegments(double start_time, double grow_time)
{
	//update the probability distributions needed
	UpdateToTime(start_time);

	const SegmentTipList & tiplist = culture->GetActiveTipList();
	auto iter = tiplist.begin();
	while (iter != tiplist.end())
	{
		GrowSegment(*iter, start_time, grow_time);
		++iter;
	}

}
//calculates the time a segment took to grow during the current timestep, this assumes the growthrate is constant within a step
double FragmentBranching::TimeOfGrowth(Segment * seg, double time_of_growth)
{
	assert(time_of_growth != 0.0);
	return time_of_growth * (seg->length()/seg->expected_length);
}

FragmentBranching * FragmentBranching::GetBrancherForSegment(Segment * seg)
{
	FragmentBranching * fb = nullptr;
	//try the caller first
	if (std::find(this->culture->m_pmat->domainptrs.begin(), this->culture->m_pmat->domainptrs.end(), seg->tip(1).pt.ndomain) != this->culture->m_pmat->domainptrs.end())
	{
		//change this if tips can grow across materials
		fb = this;
	}
	else
	{
		//find the material that contains this domain
		for (size_t j = 0; j < this->culture->m_pmat->m_pangio->m_pmat.size(); j++)
		{
			if (std::find(this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.begin(), this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.end(),
				seg->tip(1).pt.ndomain) != this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.end())
			{
				fb = this->culture->m_pmat->m_pangio->m_pmat[j]->GetCommonAngioProperties()->fbrancher;
				break;//only one angio material should contain a domain
			}
		}
	}
	assert(fb);//should only grow into valid materials
	return fb;
}

void DirectionalWeights(double da, double dw[2])
{
	double degree_isotropy = 1.0 - da;
	dw[0] = degree_isotropy;
	dw[1] = da;
}

void NoFragmentBranching::GrowSegment(Segment::TIP * tip, double starttime, double grow_time)
{
	Segment seg = culture->GrowSegment(*tip, starttime, grow_time);
	culture->AddNewSegment(seg);
}




//begin implementation of psuedo deferred branching
PsuedoDeferedFragmentBranching::PsuedoDeferedFragmentBranching(FEModel * model) : FragmentBranching(model)
{
	AddProperty(&length_to_branch_point, "length_to_branch");
	AddProperty(&time_to_emerge, "time_to_emerge");
}
void PsuedoDeferedFragmentBranching::GrowSegment(Segment::TIP * tip, double start_time, double grow_time)
{
	if (!tip->bactive)
		return;
	assert(tip->length_to_branch > 0.0);
	assert(tip->wait_time_to_branch >= 0.0);
	static int nseg_add = 0;
	Segment seg = culture->GrowSegment(*tip, start_time, grow_time);
	if (seg.length() < culture->m_pmat->m_cultureParams.min_segment_length)
	{
		printf("segment ended due to length\n");
		tip->bactive = true;
		return;
	}
	seg.expected_length = seg.length();
		
	//now calculate the length to branch from the new tip
	seg.tip(1).length_to_branch = tip->length_to_branch;
	seg.tip(0).length_to_branch = tip->length_to_branch;
	seg.tip(1).wait_time_to_branch = tip->wait_time_to_branch;
	
	culture->AddNewSegment(seg);
	auto rseg = culture->RecentSegments();
	if (rseg.size() == 1) //handle segments in non bouncing mode
	{
		rseg[0]->tip(1).length_to_branch -= seg.length();
		rseg[0]->SetTimeOfBirth(start_time);
		rseg[0]->tip(0).length_to_branch = tip->length_to_branch;
		rseg[0]->tip(1).length_to_branch = tip->length_to_branch - rseg[0]->length();
		rseg[0]->tip(0).wait_time_to_branch = tip->wait_time_to_branch;
		rseg[0]->tip(1).wait_time_to_branch = tip->wait_time_to_branch;
		tip->connected = rseg[0];
		rseg[0]->tip(0).connected = tip->parent;
		if (rseg[0]->tip(1).length_to_branch < 0.0)
		{
			double bf = (rseg[0]->length() + rseg[0]->tip(1).length_to_branch) / rseg[0]->length();
			assert(bf >= 0.0 && bf <= 1.0);
			double bt = mix(start_time, start_time + TimeOfGrowth(rseg[0], grow_time), bf);
			//add the points
			//change this if tips can grow across materials
			BranchPoint bp(bt + tip->wait_time_to_branch, bt, rseg[0], bf, rseg[0]->m_nid, GetBrancherForSegment(rseg[0]));
			branch_points.insert(bp);
			parentgen.insert(bp);
		}
	}
	else if (rseg.size() == 0)
	{
		//reactivate the tip
		tip->bactive = true;
		nseg_add++;
		//printf("segment reactivated\n");
	}
	else
	{
		for (size_t i = 0; i < rseg.size(); i++)
		{
			if (rseg.size() >(i + 1))
			{
				rseg[i]->tip(1).connected = rseg[i + 1];
				rseg[i + 1]->tip(0).connected = rseg[i];
			}
		}
		assert(seg.tip_c(1).length_to_branch > 0.0);//force the update/intialization of values
		assert(seg.tip_c(1).wait_time_to_branch >= 0.0);
		double ct = start_time;
		double l2b = seg.tip_c(1).length_to_branch;
		tip->connected = rseg[0];
		rseg[0]->tip(0).connected = tip->parent;
		//set all of the adjacenty pointers
		
		for (size_t i = 0; i < rseg.size(); i++)
		{
			rseg[i]->SetTimeOfBirth(ct);
			rseg[i]->tip(0).length_to_branch = l2b;
			l2b -= rseg[i]->length();
			rseg[i]->tip(1).length_to_branch = l2b;
			rseg[i]->tip(1).wait_time_to_branch = seg.tip_c(1).wait_time_to_branch;
			assert(rseg[i]->tip_c(0).connected);
			if (l2b <= 0.0)
			{
				double bf = (rseg[i]->length() + rseg[i]->tip_c(1).length_to_branch) / rseg[i]->length();
				assert(bf >= 0.0 && bf <= 1.0);
				double bt = mix(ct, ct + TimeOfGrowth(rseg[i], grow_time), bf);
				FragmentBranching * fb = nullptr;
				//add the points
				assert(tip->wait_time_to_branch >= 0.0);
				//note this is what GetBrancherForSegment is based on but also contains logic for BranchPoints
				if (std::find(this->culture->m_pmat->domainptrs.begin(), this->culture->m_pmat->domainptrs.end(), rseg[i]->tip_c(1).pt.ndomain) != this->culture->m_pmat->domainptrs.end())
				{
					//change this if tips can grow across materials
					BranchPoint bp(bt + tip->wait_time_to_branch, bt, rseg[i], bf, rseg[i]->m_nid, this);
					fb = this;
					branch_points.insert(bp);
					parentgen.insert(bp);
				}
				else
				{
					//find the material that contains this domain
					for (size_t j = 0; j < this->culture->m_pmat->m_pangio->m_pmat.size(); j++)
					{
						if (std::find(this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.begin(), this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.end(),
							rseg[i]->tip_c(1).pt.ndomain) != this->culture->m_pmat->m_pangio->m_pmat[j]->domainptrs.end())
						{
							fb = this->culture->m_pmat->m_pangio->m_pmat[j]->GetCommonAngioProperties()->fbrancher;
							//change this if tips can grow across materials
							BranchPoint bp(bt + tip->wait_time_to_branch, bt, rseg[i], bf, rseg[i]->m_nid, fb);
							branch_points.insert(bp);
							parentgen.insert(bp);
							break;//only one angio material should contain a domain
						}
					}
				}
				assert(fb);
				break;//only insert a single branch point
			}

			ct += TimeOfGrowth(rseg[i], grow_time);
		}
	}
}

void PsuedoDeferedFragmentBranching::GrowSegment(std::set<BranchPoint>::iterator bp,double start_time, double grow_time)
{
	//create a tip of origin
	static int mis_count = 0;
	Segment::TIP tip0 = bp->parent->tip(0);
	Segment::TIP tip1 = bp->parent->tip(1);
	vec3d t0p = bp->parent->tip_c(0).pos();
	vec3d t1p = bp->parent->tip_c(1).pos();
	vec3d pos = mix(t0p, t1p , bp->percent_of_parent);
	if (culture->m_pmat->FindGridPoint(pos, tip0.pt.ndomain, tip0.pt.elemindex, tip0.pt))
	{
		//roll for the length to branch of the new segment
		//the overwrite here is okay as it is a copy
		tip0.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
		tip0.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
		assert(tip0.wait_time_to_branch > 0.0);
		culture->BranchSegment(tip0, bp->emerge_time, start_time + grow_time - bp->emerge_time);
	}
	else if (culture->m_pmat->FindGridPoint(pos, tip1.pt.ndomain, tip1.pt.elemindex, tip1.pt))
	{
		tip1.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
		tip1.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
		assert(tip1.wait_time_to_branch > 0.0);
		culture->BranchSegment(tip1, bp->emerge_time, start_time + grow_time - bp->emerge_time);
	}
	else
	{
		//this is when the segment is in another element ie the segment is growing through at least 2 elements
		//this may indicate a problem in getting consistent results
		if (culture->m_pmat->FindGridPoint(pos, tip0.pt))
		{
			tip0.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
			tip0.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
			assert(tip0.wait_time_to_branch > 0.0);
			culture->BranchSegment(tip0, bp->emerge_time, start_time + grow_time - bp->emerge_time);
		}
		else
		{
			//need to try other materials
			bool found = false;
			for (size_t i = 0; i < culture->m_pmat->m_pangio->m_pmat.size(); i++)
			{
				//this is when the segment is in another element ie the segment is growing through at least 2 elements
				//this may indicate a problem in getting consistent results
				//don't retest the material that has already been done
				if ((culture->m_pmat != culture->m_pmat->m_pangio->m_pmat[i]) && culture->m_pmat->m_pangio->m_pmat[i]->FindGridPoint(pos, tip0.pt))
				{
					tip0.length_to_branch = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
					tip0.wait_time_to_branch = time_to_emerge->NextValue(culture->m_pmat->m_pangio->rengine);
					assert(tip0.wait_time_to_branch > 0.0);
					culture->m_pmat->m_pangio->m_pmat[i]->m_cult->BranchSegment(tip0, bp->emerge_time, start_time + grow_time - bp->emerge_time);
					found = true;
					break;
				}
			}
			//failed to place the tip in another element
			//assert(found);
		}
	}
	ProcessNewSegments(bp->emerge_time,grow_time);

}

void PsuedoDeferedFragmentBranching::ProcessNewSegments(double start_time, double grow_time)
{
	auto rseg = culture->RecentSegments();

	//consider handling when no segments are added
	if (rseg.size() == 0)
	{
		printf("no segments added\n");
	}
	//the easy case just a single segment has been added
	else
	{
		double ct = start_time;
		double lr = rseg[0]->tip_c(0).length_to_branch;
		assert(lr > 0.0);
		for (size_t i = 0; i < rseg.size(); i++)
		{
			rseg[i]->SetTimeOfBirth(ct);
			ct += TimeOfGrowth(rseg[i],grow_time);
			if (i)
			{
				//set the connections
				rseg[i - 1]->tip(1).connected = rseg[i];
				rseg[i]->tip(0).connected = rseg[i - 1];
			}
		}
		for (size_t i = 0; i < rseg.size(); i++)
		{
			rseg[i]->tip(0).length_to_branch = lr;
			lr -= rseg[i]->length();
			rseg[i]->tip(1).length_to_branch = lr;
			if (lr < 0.0)
			{
				//add a Branch point and break out
				double len = rseg[i]->length();

				double bpct = (len + lr) / len;
				lr += length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
				rseg[i]->tip(1).length_to_branch = lr;
				//add it back to the queue
				double bt = mix(rseg[i]->GetTimeOfBirth(), rseg[i]->GetTimeOfBirth() + TimeOfGrowth(rseg[i],grow_time), bpct);
				//consider doing more to determine which material the tip is in
				BranchPoint bpt(bt + rseg[0]->tip_c(1).wait_time_to_branch, bt, rseg[i], bpct, rseg[i]->m_nid, GetBrancherForSegment(rseg[i]));
				parentgen.insert(bpt);
				//make sure adding both back was wrong
				branch_points.insert(bpt);
				break;
			}
		}
		//currently branching in the initial growth of a branch is not implemented
		assert(lr > 0.0);
	}
}
void PsuedoDeferedFragmentBranching::UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp, double start_time, double grow_time)
{
	//check that the GetBrancherForSegment work properly they may be 1 segment late on switching 
	double diff = length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
	bp->parent->tip(1).length_to_branch += diff;
	bp->parent->tip(0).length_to_branch += diff;
	if (bp->parent->tip(1).length_to_branch < 0.0)
	{
		assert(bp->parent->tip_c(0).length_to_branch >= 0.0);
		double bpct = bp->parent->tip_c(0).length_to_branch / bp->parent->length();
		//add it back to the queue
		double bt = mix(bp->epoch_time, start_time + grow_time, bpct);
		
		BranchPoint bpt(bt + bp->parent->tip_c(1).wait_time_to_branch, bt, bp->parent, bpct, bp->priority, GetBrancherForSegment(bp->parent));
		parentgen.insert(bpt);
		//make sure adding both back were wrong
		branch_points.insert(bpt);
		return;
	}
	//propogate the length to branch point 
	Segment * cur = bp->parent->tip_c(1).connected;
	double dr = bp->parent->tip_c(1).length_to_branch;
	//may need to have an additional condition when anastamosis is added or anastamosis may not connect segments
	while(cur)
	{
		double len = cur->length();
		dr = dr - len;
		cur->tip(1).length_to_branch = dr;
		if (dr < 0.0)
		{
			double bpct = (len + dr) / len;
			dr += length_to_branch_point->NextValue(culture->m_pmat->m_pangio->rengine);
			cur->tip(1).length_to_branch = dr;
			//add it back to the queue
			double bt = mix(bp->epoch_time, bp->epoch_time + TimeOfGrowth(cur,grow_time), bpct);
			
			BranchPoint bpt(bt + bp->parent->tip_c(1).wait_time_to_branch, bt, bp->parent, bpct, bp->priority, GetBrancherForSegment(bp->parent));
			parentgen.insert(bpt);
			//make sure adding both back was wrong
			branch_points.insert(bpt);
			break;//only one point should be added at a time
		}

		cur = cur->tip(1).connected;
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

void PsuedoDeferedFragmentBranching::UpdateToTime(double starttime)
{
	time_to_emerge->StepToTime(starttime);
	length_to_branch_point->StepToTime(starttime);
}