#pragma once

#include <list>
#include "BC.h"
using namespace std;

//-----------------------------------------------------------------------------
class Grid;
class Segment;
class FEAngio;

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

	// create initial segments
	void SeedFragments();

public:
	// Create a new segment at the tip of an existing segment
	Segment createNewSeg(list<Segment>::iterator it, int k);
	
	// Determine the orientation angle of a newly created segment based on the information stored in GRID
	vec3d findAngle(list<Segment>::iterator, double xpt, double ypt, double zpt);
	
	// Obtain the component of new vessel orientation determined by local collagen fiber orientation 
	vec3d findCollAngle(double xpt, double ypt, double zpt);
	
	// Find the density scale factor at a point of the grid
	// TODO: Move to Grid class
	double findDenScale(double xpt, double ypt, double zpt);
	
	// Create a new segment connecting two existing segments that are fusing through anastomosis
	Segment connectSegment(list<Segment>::iterator it, list<Segment>::iterator it2, int k, int kk);
	
	// Check a newly created segment to see if it physically intersections with any existing segments
	void CheckForIntersection(Segment &seg, list<Segment>::iterator it);
	
	// Find the coordinates at which two segments intersect
	double findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3]);
	
	// Determine if a segment encounters one of the boundary planes, find the coordinates of the intersection point
	bool intersectPlane(Segment &Seg, int n, double intersectpt[3]);
	
	// If a segment encounters one of the boundary planes, enforce the periodic boundary conditions
	Segment PeriodicBC(Segment &seg);

private:
    // Seed the initial fragments within the domain
	Segment createInitFrag();

public: // TODO: make private
	list<Segment> m_frag;		// vessel fragments

public:
    BC		bc;
	double	W[4];			// W[1] = Weight for vessel density
    int		NFRAGS;			// Number of initial microvessel fragments
    int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments

	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

private:
	FEAngio&	m_angio;
};
