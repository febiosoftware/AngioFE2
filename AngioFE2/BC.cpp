#include "stdafx.h"
#include <iostream>
#include "BC.h"
#include "Culture.h"
#include "Segment.h"
#include "FEAngio.h"
#include "angio3d.h"

//-----------------------------------------------------------------------------
BC::BC(FEAngio& angio) : m_angio(angio)
{
}                                   

//-----------------------------------------------------------------------------
BC::~BC()
{

}

//-----------------------------------------------------------------------------
// Checks if a newly-created segment violates the boundary faces of the element in which it occupies
void BC::CheckBC(Segment &seg)
{
	Grid& grid = m_angio.GetGrid();

	// get the end-points and reference element
	vec3d r0 = seg.tip(0).pos();
	vec3d r1 = seg.tip(1).pos();
	int elem_num = seg.tip(0).pt.nelem;
	assert(elem_num >= 0);

	// find the intersection with the element's boundary
	FACE_INTERSECTION ic;
	ic.nelem = elem_num;
	if (grid.FindIntersection(r0, r1, ic))
	{
		// enforce the BC
		EnforceBC(seg, ic);
	}
	else
	{
		assert(false);
	}
} 

//-----------------------------------------------------------------------------
// This enforces a boundary condition on a new vessel.
void BC::EnforceBC(Segment &seg, FACE_INTERSECTION& ic)
{
	Grid& grid = m_angio.GetGrid();
	Culture& cult = m_angio.GetCulture();

	// get the element
	Elem& elem = grid.GetElement(ic.nelem);

	// get the BC type
	unsigned int bctype = grid.GetFace(ic.nface).bc_type;
	
	// vessel stops growing
    if (bctype == BC::STOP){
        BCStop(seg, ic);
		return;
	}

	// vessl bounces off the wall
	if (bctype == BC::BOUNCY){
		BCBouncy(seg, ic);
		return;
	}

	assert(false);
/*    
    // Bouncy wall boundary type
    if (bctype == 98){
        Segment seg2;
		
		seg2 = bouncywallBC(i_point, face, seg, elem_num, k);

		elem_num = grid.findelem(seg2.tip(k).rt);
        
		if (elem_num != -1){
			seg2.tip(k).pt.nelem = elem_num;
			cult.AddSegment(seg2);}
		else{
			checkBC(seg2, k);
			elem_num = grid.findelem(seg2.tip(k).rt);
			if (elem_num != -1){
				seg2.tip(k).pt.nelem = elem_num;
				cult.AddSegment(seg2);}}

		BC_bouncy = false;
		BC_violated = false;
		
		if (seg2.mark_of_death == true)
			seg.tip(k).BC = 1;

		return;}
    
	// In-plane wall boundary type
    if (bctype == 105){
        Segment seg2;
				
		seg2 = inplanewallBC(i_point, face, seg, elem_num, k);

		elem_num = grid.findelem(seg2.tip(k).rt);
        
		if (elem_num != -1){
			seg2.tip(k).pt.nelem = elem_num;
			cult.AddSegment(seg2);}
		else{
			checkBC(seg2, k);
			elem_num = grid.findelem(seg2.tip(k).rt);
			if (elem_num != -1){
				seg2.tip(k).pt.nelem = elem_num;
				cult.AddSegment(seg2);}}

		BC_bouncy = false;
		BC_violated = false;
		
		if (seg2.mark_of_death == true)
			seg.tip(k).BC = 1;

		return;} 

	// Collagen Fiber Bouncy wall boundary type
    if (bctype == 99){
        collfibwallBC(i_point, face, seg, elem_num, k);
        BC_violated = false;
		return;}   

	// Sym plane periodic wall boundary type
    if (bctype == 112){
        Segment seg2;
				
		seg2 = symplaneperiodicwallBC(i_point, face, seg, elem_num, k);

		elem_num = grid.findelem(seg2.tip(1).rt);
        
		if (elem_num != -1){
			seg2.tip(k).pt.nelem = elem_num;
			cult.AddSegment(seg2);}
		else{
			checkBC(seg2, k);
			elem_num = grid.findelem(seg2.tip(k).rt);
			if (elem_num != -1){
				seg2.tip(k).pt.nelem = elem_num;
				cult.AddSegment(seg2);}}

		BC_bouncy = false;
		BC_violated = false;
		
		if (seg2.mark_of_death == true)
			seg.tip(k).BC = 1;

		return;} 
*/
}

//-----------------------------------------------------------------------------
// Boundary condition where vessels stops growing after hitting boundary.
void BC::BCStop(Segment &seg, FACE_INTERSECTION& ic)
{
	Segment::TIP& tip = seg.tip(1);
	assert(tip.pt.nelem == -1);

	// update position and grid point structure
	Grid& grid = m_angio.GetGrid();
	tip.pt.nelem = ic.nelem;
	vec3d p = ic.q;
	grid.natcoord(tip.pt.q, p, ic.nelem);
	tip.pt.r = ic.q;
	seg.Update();

	// turn the tip off
	tip.bactive = false;

	// mark as dead
	seg.SetFlagOn(Segment::BC_DEAD);
	tip.BC = 1;

	// add the segment
	Culture& cult = m_angio.GetCulture();
	cult.AddSegment(seg);
}

//-----------------------------------------------------------------------------
// Boundary condition where the vessel bounces off the wall.
// This effectively creates another segment by breaking the current segment
// in two at the intersection point and then creating a new segment that bounces of the wall.
void BC::BCBouncy(Segment &seg, FACE_INTERSECTION& ic)
{
	Grid& grid = m_angio.GetGrid();
	Culture& cult = m_angio.GetCulture();

	// get the original tip positions
	vec3d r0 = seg.tip(0).pos();
	vec3d r1 = seg.tip(1).pos();
	vec3d q = ic.q;	// intersection point

	// get the old length
    double old_length = seg.length();
	if (old_length == 0.0)
	{
//		assert(false);
		return;
	}

	// break the first segment
	Segment::TIP& tip = seg.tip(1);
	tip.pt.r = q;
	tip.pt.nelem = ic.nelem;
	grid.natcoord(tip.pt.q, q, ic.nelem);
	seg.Update();
	tip.bactive = false;
    seg.SetFlagOn(Segment::BC_DEAD);

	// calculate new length
	double new_length = seg.length();

	// Add this segment
	// TODO: if the new length is zero, I should not add it, but then the
	//       continuity of the vessel will be broken. Not sure yet how to handle this.
	cult.AddSegment(seg);

	// create a new segment
	Segment seg2;
	seg2.seed(seg.seed());
	seg2.vessel(seg.vessel());

	// set the starting point at the intersection point
	seg2.tip(0).pt = tip.pt;
	seg2.tip(0).bactive = false;

	// calculate the recoil vector
	vec3d t = r1 - r0; t.unit();
	vec3d n = ic.norm;

	vec3d new_vec = t - n*(2.0*(t*n));
	
	// calculate the remaining length
    double remain_length = old_length - seg.length();
	assert(remain_length > 0.0);

	seg2.tip(1).pt.r = q + new_vec*remain_length;
	seg2.Update();

	// activate the new tip
	seg2.tip(1).bactive = true;

	// pass on body force ID
	seg2.tip(1).bdyf_id = seg.tip(1).bdyf_id; 

	// Add the new segment to the culture
	cult.AddNewSegment(seg2);
}
