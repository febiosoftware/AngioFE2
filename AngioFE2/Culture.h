#pragma once

#include <list>
#include "BC.h"
using namespace std;

//-----------------------------------------------------------------------------
class Grid;
class Segment;
class FEAngio;
class SimulationTime;
class GridPoint;

//-----------------------------------------------------------------------------
typedef list<Segment>::iterator SegIter;

//-----------------------------------------------------------------------------
// The CULTURE class contains all the functions that describe how the 
// SEGMENT class is to grow and orient itself. These functions are the 
// rule-set that arrange the line segments into vascular networks, 
// mimicking angiogenesis. 
class Culture  
{
public:
	Culture(FEAngio& angio);
	virtual ~Culture();

	// initialize
	bool Init();

	// create initial segments
	void SeedFragments(SimulationTime& time);

	// Perform a growth step
	void Grow(SimulationTime& time);

	// perform a sub-growth step
	void SubGrowth(double scale);

public:
	// Determine the orientation vector of a newly created segment
	vec3d FindDirection(Segment& it, GridPoint& pt);
	
	// find the unit direction vector of the collagen
	vec3d CollagenDirection(GridPoint& pt);
	
	// Find the density-based length scale factor at a point of the grid
	double FindDensityScale(const GridPoint& pt);

	// scales the density
	double ScaleDensity(double coll_den);
	
	// Check a newly created segment to see if it physically intersections with any existing segments
	void CheckForIntersection(Segment &seg, Segment& it);
	
public:
	// Add a segment to the culture
	void AddSegment(Segment& seg);

	// return the number of segments
	int Segments() { return m_nsegs; }

	// get a segment iterator
	SegIter SegmentBegin() { return m_frag.begin(); }

	// get the end iterator
	SegIter SegmentEnd() { return m_frag.end(); }

private:
    // Seed an initial fragment within the grid
	Segment createInitFrag();

	// create a branch
	void Branch(Segment& it, SimulationTime& time);

	// fuse segments (i.e. anastomosis)
	void Fuse(SimulationTime& time);

	// Create a new segment at the tip of an existing segment
	Segment CreateNewSeg(Segment& it, int k, SimulationTime& time, bool branch = false);

	// Create a new segment connecting two existing segments that are fusing through anastomosis
	Segment ConnectSegment(Segment& it, Segment& it2, int k, int kk, SimulationTime& time);

	// Update the new vessel length 
	void UpdateNewVesselLength(SimulationTime& time);

	// Find the active tips
	void FindActiveTips();

	void CreateBranchingForce(Segment& seg);
    void check4anast(Segment& it, int k, SimulationTime& time);
    void anastomose(double dist0, double dist1, int k, Segment& it, Segment& it2, SimulationTime& time);


	void kill_dead_segs();

public: // TODO: make private
	list<SegIter> m_active_tips;		// list of active segments

private:
	int				m_nsegs;	// Counter that stores in current number of Segments within the simulation domain
	list<Segment>	m_frag;		// vessel fragments

public:
    BC		bc;
	double	W[4];			// W[1] = Weight for vessel density
    int		m_ninit_frags;	// Number of initial microvessel fragments
    int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments
	int		m_num_branches;		// Counter indicating the number of branches formed during the simulation
    
	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 	                                                            
	// TODO: Makes these user parameters
	//          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double	m_y0;   // Bottom of sigmoid curve
	double	m_a;    // Distance between top and bottom of the curve (a + y0 = top of curve)
	double	m_x0;   // Time point at which t is halfway between top & bottom
	double	m_b;    // Steepness of the curve
	double	m_d;    // Initial value of the curve (t = 0)
	double	m_vess_length;		// new segment length

    double	m_init_branch_prob;     // Probability that initial segments will branch (TODO: Make this a user parameter)

    // Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjust;                                       

	// If the shortest distance vector between a segment and any other segment is less than this value
	// TODO: I think making this a percentage of the growth length makes more sense.
    double m_anast_dist;

	double	m_branch_chance;    // Probability of forming a new branch
	bool	yes_branching;
	bool	yes_anast;

	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

private:
	FEAngio&	m_angio;
};
