#pragma once
#include "StdAfx.h"
#include "BC.h"
#include "Segment.h"
#include "CultureParameters.h"
#include "FragmentSeeder.h"
#include "FragmentBranching.h"
#include "GrowDirectionModifier.h"


//-----------------------------------------------------------------------------
class FEAngio;
class SimulationTime;
class GridPoint;
class FEAngioMaterial;
class BouncyBC;

//-----------------------------------------------------------------------------
typedef list<Segment> SegmentList;
typedef list<Segment>::iterator SegIter;
typedef list<Segment>::const_iterator ConstSegIter;
typedef list<Segment::TIP*>	SegmentTipList;
typedef list<Segment::TIP*>::iterator TipIter;
typedef list<Segment::TIP*>::const_iterator ConstTipIter;






//-----------------------------------------------------------------------------
// The CULTURE class contains all the functions that describe how the 
// SEGMENT class is to grow and orient itself. These functions are the 
// rule-set that arrange the line segments into vascular networks, 
// mimicking angiogenesis. 
class Culture  
{
public:
	Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp, FragmentBranching *fbr);
	virtual ~Culture();

	// initialize
	bool Init();

	// Reposition the vessels based on the FE solution
	void Update();

	// get the total vessel length
	double TotalVesselLength() const { return m_total_length; }
	
	// Find the density-based length scale factor at a point of the grid
	double FindDensityScale(const GridPoint& pt) const;
	
	// Add a segment to the culture (without checking BCs).
	// This assumes that the segment is valid.
	void AddSegment(Segment& seg);

	// Add a new segment to the culture.
	// This will apply BCs to the new segment and may result in 
	// the addition of several new segments. 
	void AddNewSegment(Segment& seg);

	//same as AddNewSegment but does not clear the recent segments used in bouncy boundaries
	void AddNewSegmentNoClear(Segment& seg);

	//returns the segments that were added since the last call to AddNewSegment
	const std::vector<Segment *> & RecentSegments() const { return recents; }

	// return the number of segments
	int Segments() const { return m_nsegs; }

	// get the segment list
	const SegmentList& GetSegmentList() const { return m_frag; }

	// return the active segment list
	const SegmentTipList& GetActiveTipList() const
	{
		return m_active_tips;
	}
	//returns the active tips sorted by x position, the sorting key is arbitrary but should give consistent results between runs
	//consider sorting by distance from the origin instead
	const SegmentTipList & GetActiveSortedTipList()
	{
		m_active_tips.sort([](Segment::TIP * t0, Segment::TIP * t1){return t0->pos().x < t1->pos().x; });
		return m_active_tips;
	}

	// get the total number of active tips
	int ActiveTips() const { return static_cast<int>(m_active_tips.size()); }

	void CreateBranchingForce(Segment& seg);

	// Find the active tips
	void FindActiveTips();

	// Grow a segment
	Segment GrowSegment(Segment::TIP& it,double starttime, double grow_time, bool branch = false, bool bnew_vessel = false);

	//returns the segment length given the time parameters this will be further adjusted by ECM density
	double SegmentLength(double starttime, double grow_time) const;

	// create a branch
	void BranchSegment(Segment::TIP& it, double starttime, double grow_time);


private:	
	// Create a new segment connecting two existing segments that are fusing through anastomosis
	static Segment ConnectSegment(Segment& it, Segment& it2, int k, int kk);

	


	int					m_nsegs;				// Counter that stores in current number of Segments within the simulation domain
	double				m_total_length;			// Total vascular length within the domain (sum of the length of all Segments) (in um)
	list<Segment>		m_frag;					// vessel fragments
	list<Segment::TIP*> m_active_tips;			// list of active tips
	CultureParameters * m_cultParams;
	FEAngio&	m_angio;
	
	std::vector<Segment *> recents;//used to hold the segments added by the most recent call to AddNewSegment these segments will be ordered first to last

public:
	FragmentBranching * fbrancher = nullptr;
	double	m_vess_length;	// new segment length
	int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments
	int		m_num_branches;		// Counter indicating the number of branches formed during the simulation
	FEAngioMaterial * m_pmat = nullptr;
	
	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

	

	friend class BouncyBC;
	friend class BC;
	friend class StopBC;
	friend class MBC;
	friend class PassThroughMBC;
};
