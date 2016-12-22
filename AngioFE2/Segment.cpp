#include "StdAfx.h"
#include "Segment.h"
#include <FECore/vec3d.h>

//-----------------------------------------------------------------------------
Segment::TIP::TIP()
{
	bactive = false;
	BC = 0;
	nseed = -1;
	nvessel = -1;
}

//-----------------------------------------------------------------------------
Segment::Segment()
{
	// initialize data
	m_length = 0;                                        
	m_TofBirth = 0;
	
	// initialize IDs
	m_nid = 0;
	m_nseed   = -1;
    m_nvessel = -1;
    
	// initialize status flags
	m_nflag = 0;

	tip(0).parent = this;
	tip(1).parent = this;
}
Segment::Segment(const Segment &obj)
{
	//can i just do what ever would have been done before for the rest of the parameters
	memcpy(this, &obj, sizeof(Segment));
	//update the parent pointers
	tip(0).parent = this;
	tip(1).parent = this;
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
	m_uvect = m_tip[1].pos() - m_tip[0].pos();
	m_length = m_uvect.unit();

	m_tip[0].u = -m_uvect;
	m_tip[1].u =  m_uvect;
}
