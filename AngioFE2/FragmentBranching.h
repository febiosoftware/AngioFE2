#pragma once

#include "StdAfx.h"
#include "Segment.h"
#include "FECore/FEMaterial.h"
#include "FEProbabilityDistribution.h"

class Culture;


//this class encapsulates the branching behavior for a given culture(instance)
//the static methods/members synchronize the instances where needed, these should only be called once from the model
class FragmentBranching : public FEMaterial
{
public:
	class BranchPoint
	{
	public:
		BranchPoint(double emt, double ept, Segment * p, double pctp, int prior, FragmentBranching * fb) :emerge_time(emt),
			epoch_time(ept), parent(p), percent_of_parent(pctp), priority(prior), brancher(fb)
		{
			assert(emerge_time >= -1.0);
			assert(epoch_time >= -1.0);
			assert(emerge_time >= epoch_time);
			assert(parent != nullptr);
			assert(brancher != nullptr);
			//check that the parent has the l2b set
			assert(parent->tip(1).length_to_branch != 0.0);
			assert(parent->tip(0).length_to_branch != 0.0);
		}

		~BranchPoint(){}

		double emerge_time;//time at which this branch begins to grow
		double epoch_time;//time when this branch point was created
		Segment * parent;
		double percent_of_parent;//0-1.0 determines where on the parent the branch will sprout from 0 is 100% contribution from tip(0)
		int priority;//it there is a tie in time this will break it consistently needs to be athe same for the brnachpoints between runs with equivalent paramters and unique among branch points
		FragmentBranching * brancher;//used to get the rng needed for this segment

		//include utility if the other relational operators are needed
		//should allos the set to be iterated over from low to high times
		bool operator< (const BranchPoint& rhs) const
		{
			if (emerge_time < rhs.emerge_time)
				return true;
			else if (emerge_time > rhs.emerge_time)
				return false;
			else
				return priority < rhs.priority;
		}
		bool operator== (const BranchPoint & rhs) const
		{
			return emerge_time == rhs.emerge_time && parent == rhs.parent && percent_of_parent == rhs.percent_of_parent && priority == rhs.priority;
		}

	};
	//allows Branch points to be sorted w/i a set by epoch time
	struct BranchPointEpochCompare
	{
		bool operator()(const BranchPoint & lhs, const BranchPoint & rhs) const
		{
			if (lhs.epoch_time < rhs.epoch_time)
				return true;
			else if (lhs.epoch_time > rhs.epoch_time)
				return false;
			else
				return lhs.priority < rhs.priority;
		}
	};
	//allows Branch points to be sorted w/i set sorted by the lowest of the times
	struct BranchPointTimeFloorCompare
	{
		bool operator()(const BranchPoint & lhs, const BranchPoint & rhs) const
		{
			return (std::min(lhs.epoch_time, lhs.emerge_time) < std::min(rhs.epoch_time, rhs.emerge_time));
		}
	};


	FragmentBranching(FEModel * model):FEMaterial(model)
	{
		fragment_branchers.push_back(this);
	}
	
	virtual ~FragmentBranching()
	{
		assert(std::find(fragment_branchers.begin(), fragment_branchers.end(), this) != fragment_branchers.end());
		fragment_branchers.erase(std::find(fragment_branchers.begin(), fragment_branchers.end(), this));
	}
	//Grows a segment from a branchpoint
	virtual void GrowSegment(std::set<BranchPoint>::iterator bp) = 0;

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);

	void GrowSegments();
	//Grows a segment from an active tip
	//this function should be time independent
	//thsi grows a segment from an active tip and sets the BranchPoints needed to generate any needed rng
	virtual void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) = 0;

	virtual void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) = 0;

	//modify the length to brnach if needed and set the time of birth of the segment
	virtual void PostProcess(Segment & seg) = 0;

	virtual double GetLengthToBranch() = 0;
	virtual double GetTimeToEmerge(){ return 0.0; }

	double TimeOfGrowth(Segment * seg);

	FragmentBranching * GetBrancherForSegment(Segment * seg);

	//Grow is a synchronized grow operation for all FragmentBranchers
	static void Grow();

	//call once to have this class record a timeline of the branch points it uses
	static void OutputTimeline(){ create_timeline = true; }

	friend class Fileout;
protected:
	Culture * culture;
	double current_time;
	

	static std::vector<FragmentBranching *> fragment_branchers;//consider making this private
	static std::set<BranchPoint> branch_points;//used for creating the branches
	static std::set<BranchPoint, BranchPointEpochCompare> parentgen;//used to generate the next branchpoint for the parent segment
	static bool create_timeline;
	static std::multiset<BranchPoint, FragmentBranching::BranchPointTimeFloorCompare> timeline;
};

//fragments determine the branch points as they grow but there is a time delay to when they start growing
class PsuedoDeferedFragmentBranching :public FragmentBranching
{
public:
	PsuedoDeferedFragmentBranching(FEModel * model);
	void GrowSegment(std::set<BranchPoint>::iterator bp) override;
	void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) override;
	void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) override;
	void PostProcess(Segment & seg) override;
	double GetLengthToBranch() override;
	double GetTimeToEmerge() override;
	void ProcessNewSegments(double start_time);
private:
	FEPropertyT<FEProbabilityDistribution> length_to_branch_point;
	FEPropertyT<FEProbabilityDistribution> time_to_emerge;
};
//for testing and replacing existing no branch functionality
class NoFragmentBranching : public FragmentBranching
{
public:
	NoFragmentBranching(FEModel * model) : FragmentBranching(model){}
	void GrowSegment(std::set<BranchPoint>::iterator bp) override{}
	void GrowSegment(Segment::TIP * tip, double starttime, double grow_time) override;
	void UpdateSegmentBranchDistance(std::set<BranchPoint>::iterator bp) override{}
	void PostProcess(Segment & seg) override {}
	double GetLengthToBranch() override{ return 1000000.0; }
};