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

public:
    Segment();
	virtual ~Segment();
    
    void findlength();
    void findunit();

public:
	vec3d	rt[2];		// current position of segment tips
    vec3d	uvect;		// unit direction vector

public:
	int tip[2];                                                 // SEGMENT.tip - Array indicating that the segment's nodes are active
    double length;                                              // SEGMENT.length - Length of the segment (in um), Euclidean distance between the nodes
	
	int label;                                                  // SEGMENT.label - Label that indicates which initial fragment the segment orginated from
	int vessel;                                                 // SEGMENT.vessel - Label that indicates which vessel the segment belongs to
	int seg_num;
	
	int BCdead;                                                 // SEGMENT.BCdead - Boolean flag that indicates that the segment has encountered a boundary condition
	double TofBirth;                                            // SEGMENT.TofBirth - Time point at which the segment was created
	double Recent_branch;                                       // SEGMENT.Recent_branch - Indicates that the segment was recently involved in the formation of a branch
	bool init_branch;                                           // SEGMENT.init_branch - Boolean flag that indicates whether or not the initial fragment should form a branch
	                                                            //                       at t = 0
	int m_sprout;                                                 // SEGMENT.sprout - Marker that indicates which end of the parent vessel the segment sprouted from.  1 for the +1 end, 
	                                                            //                  -1 for the -1 end, 9 for an initially seeded fragment,
	                                                            //                  0 for a frament not touched by the growth routine for some reason.
    int anast;
    
	int tip_elem[2];                                               // SEGMENT.seg_elem - Indicates which element the segment is currently occupying

	bool elem_tagged;

	int bdyf_id[2];

	bool mark_of_death;
	int death_label;

	int tip_BC[2];

	int seg_conn[2][2];
};
