#include "stdafx.h"
#include "Segment.h"
#include <FECore/vec3d.h>
#include "angio3d.h"

//-----------------------------------------------------------------------------
Segment::TIP::TIP()
{
	rt = vec3d(0,0,0);
	bactive = false;
	bdyf_id = -1;
	BC = 0;
}

//-----------------------------------------------------------------------------
Segment::Segment()
{
	// initialize data
	m_length = 0;                                        
	m_TofBirth = 0;
	
	// initialize IDs
	m_nid = 0;
	m_nseed = 0;
    m_nvessel = 0;
    
	// initialize status flags
	m_nflag = 0;

	// TODO: remove this
	mark_of_death = false;
	death_label = 0;
}

//-----------------------------------------------------------------------------
Segment::~Segment()
{

}

//-----------------------------------------------------------------------------
// Update the segment data (i.e. length and unit vector).
// Call this each time the position of one of the nodes has changed.
void Segment::Update()
{	
	m_length = (m_tip[1].rt - m_tip[0].rt).norm();
    
	m_uvect = m_tip[1].rt - m_tip[0].rt;
	m_uvect.unit();
}
