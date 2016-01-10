#include "stdafx.h"
#include "Segment.h"
#include <FECore/vec3d.h>
#include "angio3d.h"

//-----------------------------------------------------------------------------
Segment::TIP::TIP()
{
	rt = vec3d(0,0,0);
	active = 0;
	elem = -1;
	bdyf_id = -1;
	BC = 0;
}

//-----------------------------------------------------------------------------
Segment::Segment()
{
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

	mark_of_death = false;
	death_label = 0;
	
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
	double new_length = (m_tip[1].rt - m_tip[0].rt).norm();
    
    if (length < 0.0)
        length = -new_length;
	else
		length = new_length;
        
	uvect = m_tip[1].rt - m_tip[0].rt;
	if (length != 0)
		uvect = uvect/uvect.norm();
}


//-----------------------------------------------------------------------------
// Updates the unit vector
void Segment::findunit()
{	
	uvect = m_tip[1].rt - m_tip[0].rt;
	double l = uvect.norm();
	if (l != 0) uvect /= l;
}
