#pragma once
#include <FECore/vec3d.h>

//-----------------------------------------------------------------------------
// Microvessels are represent by a collection of line segments. 
// Growth is represented by the addition of new segments onto the 
// active tips of exisiting segments. Within the simulation, these 
// line segments are found in the SEGMENT class.
class Segment  
{
public:
	enum {
		SPROUT_UNKNOWN,	// unknown sprout type
		SPROUT_INIT,	// initial sprout
		SPROUT_POS,		// fragment sprouted from +1 end
		SPROUT_NEG		// fragment sprouted from -1 end
	};

	// struct defining a tip of the segment
	class TIP
	{
	public:
		vec3d	rt;			// current position of segment tip
		int		active;		// 1 if active otherwise 0 (I think)
		int		elem;		// the element ID the tip currently in
		int		bdyf_id;	// ID of the body force
		int		BC;			// something to do with body forces?

	public:
		TIP();
	};

public:
    Segment();
	virtual ~Segment();
    
    void findlength();
    void findunit();

public:
	TIP		m_tip[2];		// the two end tips

    double length;             // Length of the segment.
    vec3d	uvect;		// unit direction vector
	
	int label;				   // Label that indicates which initial fragment the segment orginated from
	int vessel;                // Label that indicates which vessel the segment belongs to
	int seg_num;
	
	int BCdead;                // Boolean flag that indicates that the segment has encountered a boundary condition
	double TofBirth;           // Time point at which the segment was created
	double Recent_branch;      // Indicates that the segment was recently involved in the formation of a branch
	bool init_branch;          // Boolean flag that indicates whether or not the initial fragment should form a branch at t = 0
	int m_sprout;              // Marker that indicates which end of the parent vessel the segment sprouted from.  1 for the +1 end, 
	                           //                  -1 for the -1 end, 9 for an initially seeded fragment,
	                           //                  0 for a frament not touched by the growth routine for some reason.
    int anast;				   // flag indicating the segment was merged during anastimoses.
    

	int seg_conn[2][2];			// TODO: What is this?


	// TODO: I don't think there is any reason why segments should be killed and this
	//       is probably a way to clean up after bugs. Better solution: fix bugs
	bool mark_of_death;			// this flag marks the segment as dead and should be removed.
	int death_label;			// cause of death
};
