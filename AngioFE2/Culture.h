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

public:
	// create a branch
	void Branch(list<Segment>::iterator it, SimulationTime& time);

	// fuse segments (i.e. anastomosis)
	void Fuse(SimulationTime& time);

	// Create a new segment at the tip of an existing segment
	Segment createNewSeg(list<Segment>::iterator it, int k, SimulationTime& time);

public:
	// Determine the orientation angle of a newly created segment based on the information stored in GRID
	vec3d findAngle(list<Segment>::iterator, vec3d& r);
	
	// Obtain the component of new vessel orientation determined by local collagen fiber orientation 
	vec3d findCollAngle(vec3d& pt);

	// find the unit direction vector of the collagen
	vec3d CollagenDirection(GridPoint& pt);
	
	// Find the density scale factor at a point of the grid
	// TODO: Move to Grid class
	double findDenScale(vec3d& pt);
	
	// Create a new segment connecting two existing segments that are fusing through anastomosis
	Segment connectSegment(list<Segment>::iterator it, list<Segment>::iterator it2, int k, int kk, SimulationTime& time);
	
	// Check a newly created segment to see if it physically intersections with any existing segments
	void CheckForIntersection(Segment &seg, list<Segment>::iterator it);
	
	// Find the coordinates at which two segments intersect
	double findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3]);
	
	// Determine if a segment encounters one of the boundary planes, find the coordinates of the intersection point
	bool intersectPlane(Segment &Seg, int n, double intersectpt[3]);
	
	// If a segment encounters one of the boundary planes, enforce the periodic boundary conditions
	Segment PeriodicBC(Segment &seg);

private:
    // Seed an initial fragment within the grid
	Segment createInitFrag();

	// Update the new vessel length 
	void UpdateNewVesselLength(SimulationTime& time);

	// Find the active tips
	void FindActiveTips();

	void create_branching_force(Segment& seg);
    void check4anast(list<Segment>::iterator it, int k, SimulationTime& time);
    void anastomose(double dist0, double dist1, int k, list<Segment>::iterator it, list<Segment>::iterator it2, SimulationTime& time);

	void kill_dead_segs();

public: // TODO: make private
	list<Segment> m_frag;		// vessel fragments

	list<list<Segment>::iterator > m_active_tips;		// list of active segments

public:
    BC		bc;
	double	W[4];			// W[1] = Weight for vessel density
    int		m_ninit_frags;	// Number of initial microvessel fragments
    int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments
	int		m_nsegs;		// Counter that stores in current number of Segments within the simulation domain
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

	bool	m_branch;			// Boolean flag that indicates to the model that the Segment being created is the result of a new branch
	double	m_branch_chance;    // Probability of forming a new branch
	bool	yes_branching;
	bool	yes_anast;

	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

private:
	FEAngio&	m_angio;
};
