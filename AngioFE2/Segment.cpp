#include "stdafx.h"
#include "Segment.h"
#include <FECore/vec3d.h>
#include "angio3d.h"

//-----------------------------------------------------------------------------
Segment::Segment()
{
	rt[0] = vec3d(0,0,0);
	rt[1] = vec3d(0,0,0);

	// initialize the activity of the tips to 0 (i.e. inactive)
	tip[0]= 0;
	tip[1]= 0;

	// initialize length
	length = 0;                                        
	
	label = 0;                                                  // label - Initialize label to 0
    vessel = 0;                                                 // vessel - Initialize vessel to 0
    
    BCdead = 0;                                                 // BCdead - Set boundary condition indicator to 'false'
	TofBirth = 0;                                               // TofBirth - Initialize time of birth to 0
	Recent_branch = 0;                                          // Recent_branch - Initialize branching indicator to 0
	init_branch = false;                                        // init)branch - Set initial branching flag to 'false'
	
	m_sprout = SPROUT_UNKNOWN;

    anast = 0;

	elem_tagged = false;

	bdyf_id[0] = -1;
	bdyf_id[1] = -1;
	
	mark_of_death = false;
	death_label = 0;
	
	tip_BC[0] = 0;
	tip_BC[1] = 0;

	tip_elem[0] = -1;
	tip_elem[1] = -1;

	seg_num = 0;

	seg_conn[0][0] = 0;
	seg_conn[0][1] = 0;
	seg_conn[1][0] = 0;
	seg_conn[1][1] = 0;
}

//-----------------------------------------------------------------------------
Segment::~Segment()
{

}

//-----------------------------------------------------------------------------
// Calculates the length of the segment
// TODO: This is also calculates the unit vector. Remove this.
void Segment::findlength()
{	
	double new_length = (rt[1] - rt[0]).norm();
    
    if (length < 0.0)
        length = -new_length;
	else
		length = new_length;
        
	uvect = rt[1] - rt[0];
	if (length != 0)
		uvect = uvect/uvect.norm();
}


//-----------------------------------------------------------------------------
// Updates the unit vector
void Segment::findunit()
{	
	uvect = rt[1] - rt[0];
	double l = uvect.norm();
	if (l != 0) uvect /= l;
}
