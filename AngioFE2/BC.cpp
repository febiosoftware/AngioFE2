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
	BC_violated = false;
	BC_bouncy = false;
}                                   

//-----------------------------------------------------------------------------
BC::~BC()
{

}

//-----------------------------------------------------------------------------
// Checks if a newly-created segment violates the boundary faces of the element in which it occupies
void BC::checkBC(Segment &seg, int k)
{
	Grid& grid = m_angio.GetGrid();

	// get the end-points and reference element
	vec3d r0 = (k==0? seg.tip(1).rt : seg.tip(0).rt);
	vec3d r1 = (k==0? seg.tip(0).rt : seg.tip(1).rt);
	int elem_num = (k==0? seg.tip(1).pt.nelem : seg.tip(0).pt.nelem);
	assert(elem_num >= 0);

	// find the intersection with the element's boundary
	FACE_INTERSECTION ic;
	if (grid.FindIntersection(r0, r1, elem_num, ic))
	{
		Elem& elem = grid.GetElement(elem_num);
		assert(ic.nface != -1);
		assert(elem.m_nbr[ic.nface] == -1);

		// enforce the BC
		enforceBC(seg, k, ic);
	}
	else
	{
		assert(false);
	}
} 

//-----------------------------------------------------------------------------
// This enforces a boundary condition on a new vessel.
void BC::enforceBC(Segment &seg, int k, FACE_INTERSECTION& ic)
{
	Grid& grid = m_angio.GetGrid();
	Culture& cult = m_angio.GetCulture();

	// get the element
	Elem& elem = grid.GetElement(ic.nelem);

	// get the BC type
	unsigned int bctype = elem.GetFace(ic.nface)->bc_type;
	
	// vessel stops growing
    if (bctype == BC::STOP){
        BCStop(seg, k, ic);
		return;
	}

	// vessl bounces off the wall
	if (bctype == BC::BOUNCY){
		BCBouncy(seg, k, ic);
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
void BC::BCStop(Segment &seg, int k, FACE_INTERSECTION& ic)
{
	Segment::TIP& tip = seg.tip(k);

	// update position and grid point structure
	Grid& grid = m_angio.GetGrid();
	tip.pt.nelem = ic.nelem;
	vec3d p = ic.q;
	grid.natcoord(tip.pt.q, p, ic.nelem);
	seg.tip(k).rt = ic.q;
	seg.Update();

	// turn the tip off
	seg.tip(k).bactive = false;

	// mark as dead
	seg.SetFlagOn(Segment::BC_DEAD);
	seg.tip(k).BC = 1;

	// add the segment
	Culture& cult = m_angio.GetCulture();
	cult.AddSegment(seg);
}

//-----------------------------------------------------------------------------
// Boundary condition where the vessel bounces off the wall.
// This effectively creates another segment by breaking the current segment
// in two at the intersection point and then creating a new segment that bounces of the wall.
void BC::BCBouncy(Segment &seg, int k, FACE_INTERSECTION& ic)
{
	Grid& grid = m_angio.GetGrid();
	Culture& cult = m_angio.GetCulture();

	// get the original tip positions
	vec3d r0 = seg.tip(0).rt;
	vec3d r1 = seg.tip(1).rt;
	vec3d q = ic.q;	// intersection point

	// get the old length
    double old_length = seg.length();
	if (old_length == 0.0)
	{
//		assert(false);
		return;
	}

	// break the first segment
	Segment::TIP& tip = seg.tip(k);
	tip.rt = q;
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

	// set the starting point at the intersection point
	if (k==1)
	{
		seg2.tip(0).rt = q;
		seg2.tip(0).pt = tip.pt;
		seg2.tip(0).bactive = false;
	}
	else
	{
		seg2.tip(1).rt = q;
		seg2.tip(1).pt = tip.pt;
		seg2.tip(1).bactive = false;
	}

	// calculate the recoil vector
	vec3d t = r1 - r0; t.unit();
	vec3d n = ic.norm;

	vec3d new_vec = t - n*(2.0*(t*n));
	
	// calculate the remaining length
    double remain_length = old_length - seg.length();
	assert(remain_length > 0.0);

	if (k == 1)
		seg2.tip(1).rt = q + new_vec*remain_length;
    else if (k == 0)
        seg2.tip(0).rt = q - new_vec*remain_length;

	// activate the new tip
	seg2.tip(k).bactive = true;

	// pass on body force ID
	seg2.tip(k).bdyf_id = seg.tip(k).bdyf_id; 

	// Add the new segment to the culture
	cult.AddNewSegment(seg2, k);
}

///////////////////////////////////////////////////////////////////////
// find_intersect
///////////////////////////////////////////////////////////////////////

vec3d BC::find_intersect(Elem &elem, int &face, Segment &seg)
{
	vec3d inter;
	vec3d A, B;
	vec3d R;
    mat3d Hmat;
    vec3d E;
    vec3d dE;
    vec3d newE;

	double lam = 0.5;
	double e1 = 0.;
	double e2 = 0.;

	Node X1;
	Node X2;
	Node X3;
	Node X4;

	// Check guess first
	if (face == 1){
		X1 = (*elem.n1);
		X2 = (*elem.n2);
		X3 = (*elem.n5);
		X4 = (*elem.n6);}

	if (face == 2){
		X1 = (*elem.n2);
		X2 = (*elem.n4);
		X3 = (*elem.n6);
		X4 = (*elem.n8);}

	if (face == 3){
		X1 = (*elem.n3);
		X2 = (*elem.n4);
		X3 = (*elem.n7);
		X4 = (*elem.n8);}

	if (face == 4){
		X1 = (*elem.n1);
		X2 = (*elem.n3);
		X3 = (*elem.n5);
		X4 = (*elem.n7);}

	if (face == 5){
		X1 = (*elem.n5);
		X2 = (*elem.n6);
		X3 = (*elem.n7);
		X4 = (*elem.n8);}

	if (face == 6){
		X1 = (*elem.n1);
		X2 = (*elem.n2);
		X3 = (*elem.n3);
		X4 = (*elem.n4);}

	if (seg.tip(0).bactive){
		A = seg.tip(1).rt;
		B = seg.tip(0).rt;
	}

	if (seg.tip(1).bactive){
		A = seg.tip(0).rt;
		B = seg.tip(1).rt;
	}

	newton_find_intersect(lam, e1, e2, A, B, X1, X2, X3, X4);

	Grid& grid = m_angio.GetGrid();

	if ((lam > 0) && (lam <= 1.0) && (fabs(e1) <= 1.0) && (fabs(e2) <= 1.0)){					
		inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
		inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
		inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;
				
		int elem_num = elem.elem_num;
		vec3d q(0,0,0);

		if (seg.tip(1).bactive)
			seg.tip(1).pt.nelem = elem_num;
		else if (seg.tip(0).bactive)
			seg.tip(0).pt.nelem = elem_num;
		
		grid.natcoord(q, inter, elem_num);

		if (q.x >=  1.0) q.x =  0.99;
		if (q.x <= -1.0) q.x = -0.99;
		if (q.y >=  1.0) q.y =  0.99;
		if (q.y <= -1.0) q.y = -0.99;
		if (q.z >=  1.0) q.z =  0.99;
		if (q.z <= -1.0) q.z = -0.99;
				
		grid.nattoglobal(inter, q, elem_num);

		return inter;
	}
	
	
	// If that doesn't work, check all other boundary faces
	int BC_face[6] = {0};
	if (elem.GetFace(0)) BC_face[0] = 1;
	if (elem.GetFace(1)) BC_face[1] = 1;
	if (elem.GetFace(2)) BC_face[2] = 1;
	if (elem.GetFace(3)) BC_face[3] = 1;
	if (elem.GetFace(4)) BC_face[4] = 1;
	if (elem.GetFace(5)) BC_face[5] = 1;

	for (int i = 0; i < 6; ++i)
	{
		if (BC_face[i] == 1)
		{
			face = i + 1;
		
			if (face == 1){
				X1 = (*elem.n1);
				X2 = (*elem.n2);
				X3 = (*elem.n5);
				X4 = (*elem.n6);}

			if (face == 2){
				X1 = (*elem.n2);
				X2 = (*elem.n4);
				X3 = (*elem.n6);
				X4 = (*elem.n8);}

			if (face == 3){
				X1 = (*elem.n3);
				X2 = (*elem.n4);
				X3 = (*elem.n7);
				X4 = (*elem.n8);}

			if (face == 4){
				X1 = (*elem.n1);
				X2 = (*elem.n3);
				X3 = (*elem.n5);
				X4 = (*elem.n7);}

			if (face == 5){
				X1 = (*elem.n5);
				X2 = (*elem.n6);
				X3 = (*elem.n7);
				X4 = (*elem.n8);}

			if (face == 6){
				X1 = (*elem.n1);
				X2 = (*elem.n2);
				X3 = (*elem.n3);
				X4 = (*elem.n4);}

			if (seg.tip(0).bactive){
				A = seg.tip(1).rt;
				B = seg.tip(0).rt;
			}

			if (seg.tip(1).bactive){
				A = seg.tip(0).rt;
				B = seg.tip(1).rt;
			}

			newton_find_intersect(lam, e1, e2, A, B, X1, X2, X3, X4);

	 		if ((lam > 0) && (lam <= 1.0) && (fabs(e1) <= 1.0) && (fabs(e2) <= 1.0)){					
				inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
				inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
				inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;
				
				int elem_num = elem.elem_num;
				vec3d q(0,0,0);

				if (seg.tip(1).bactive)
					seg.tip(1).pt.nelem = elem_num;
				else if (seg.tip(0).bactive)
					seg.tip(0).pt.nelem = elem_num;

				grid.natcoord(q, inter, elem_num);
				if (q.x >=  1.0) q.x =  0.99;
				if (q.x <= -1.0) q.x = -0.99;
				if (q.y >=  1.0) q.y =  0.99;
				if (q.y <= -1.0) q.y = -0.99;
				if (q.z >=  1.0) q.z =  0.99;
				if (q.z <= -1.0) q.z = -0.99;
				grid.nattoglobal(inter, q, elem_num);

				return inter;}
			else if ((fabs(e1) < 2.0) && (fabs(e2) < 2.0)){
				bool intersect_neighbor = false;
				intersect_neighbor = search_neighbors_4_intersect(elem, face, lam, e1, e2, A, B, inter);

				if (intersect_neighbor == true){
					if (seg.tip(1).bactive)
						seg.tip(1).pt.nelem = elem.elem_num;
					else if (seg.tip(0).bactive)
						seg.tip(0).pt.nelem = elem.elem_num;
					
					return inter;}}		
		}

	}

	face = -1;
	seg.death_label = 2;
	seg.mark_of_death = true;
	return inter;
}



///////////////////////////////////////////////////////////////////////
// search_neighbors_4_intersect
///////////////////////////////////////////////////////////////////////

bool BC::search_neighbors_4_intersect(Elem &elem, int face, double &lam, double &e1, double &e2, vec3d &A, vec3d &B, vec3d &inter) 
{
	Grid& grid = m_angio.GetGrid();
	int elem_num = elem.elem_num;

	Elem neighbor;
	int neighbor_num = 0;
	int righty = -1; int lefty = -1; int toppy = -1; int bottomy = -1; int fronty = -1; int backy = -1;

	double new_lam = 0.5;
	double new_e1 = 0.;
	double new_e2 = 0.;

	Node X1;
	Node X2;
	Node X3;
	Node X4;

	// Face 1
	if (face == 1){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Toppy's Lefty)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = grid.elem_find_neighbor(toppy, 4);}
			if (i == 2){		
				// Check Neighbor Face 2 (Toppy)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = toppy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Toppy's Righty)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = grid.elem_find_neighbor(toppy, 2);}
			if (i == 4){		
				// Check Neighbor Face 4 (Righty)
				righty = grid.elem_find_neighbor(elem_num, 2);
				neighbor_num = righty;}
			if (i == 5){		
				// Check Neighbor Face 5 (Bottmy's Righty)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = grid.elem_find_neighbor(bottomy, 2);}
			if (i == 6){		
				// Check Neighbor Face 6 (Bottomy)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = bottomy;}
			if (i == 7){		
				// Check Neighbor Face 7 (Bottmy's Lefty)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = grid.elem_find_neighbor(bottomy, 4);}
			if (i == 8){		
				// Check Neighbor Face 8 (Lefty)
				lefty = grid.elem_find_neighbor(elem_num, 4);
				neighbor_num = lefty;}

			
			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n1);
				X2 = (*neighbor.n2);
				X3 = (*neighbor.n5);
				X4 = (*neighbor.n6);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);

					grid.natcoord(q, inter, neighbor_num);

					if (q.x >= 1.0)
						q.x = 0.99;
					if (q.x <= -1.0)
						q.x = -0.99;
					if (q.y >= 1.0)
						q.y = 0.99;
					if (q.y <= -1.0)
						q.y = -0.99;
					if (q.z >= 1.0)
						q.z = 0.99;
					if (q.z <= -1.0)
						q.z = -0.99;
				
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);
					return true;}}}}	
	
	// Face 2
	if (face == 2){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Fronty's Toppy)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 5);}
			if (i == 2){		
				// Check Neighbor Face 2 (Toppy)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = toppy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Backy's Toppy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 5);}
			if (i == 4){		
				// Check Neighbor Face 4 (Backy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = backy;}
			if (i == 5){		
				// Check Neighbor Face 5 (Backy's Bottomy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 6);}
			if (i == 6){		
				// Check Neighbor Face 6 (Bottomy)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = bottomy;}	
			if (i == 7){		
				// Check Neighbor Face 7 (Fronty's Bottomy)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 6);}	
			if (i == 8){		
				// Check Neighbor Face 8 (Fronty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = fronty;}
			

			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n2);
				X2 = (*neighbor.n4);
				X3 = (*neighbor.n6);
				X4 = (*neighbor.n8);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);
					grid.natcoord(q, inter, elem_num);
					if (q.x >=  1.0) q.x =  0.99;
					if (q.x <= -1.0) q.x = -0.99;
					if (q.y >=  1.0) q.y =  0.99;
					if (q.y <= -1.0) q.y = -0.99;
					if (q.z >=  1.0) q.z =  0.99;
					if (q.z <= -1.0) q.z = -0.99;
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);

					return true;}}}}	

	// Face 3
	if (face == 3){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Toppy's Lefty)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = grid.elem_find_neighbor(toppy, 4);}
			if (i == 2){		
				// Check Neighbor Face 2 (Toppy)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = toppy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Toppy's Righty)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = grid.elem_find_neighbor(toppy, 2);}
			if (i == 4){		
				// Check Neighbor Face 4 (Righty)
				righty = grid.elem_find_neighbor(elem_num, 2);
				neighbor_num = righty;}
			if (i == 5){		
				// Check Neighbor Face 5 (Bottmy's Righty)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = grid.elem_find_neighbor(bottomy, 2);}
			if (i == 6){		
				// Check Neighbor Face 6 (Bottomy)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = bottomy;}
			if (i == 7){		
				// Check Neighbor Face 7 (Bottmy's Lefty)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = grid.elem_find_neighbor(bottomy, 4);}
			if (i == 8){		
				// Check Neighbor Face 8 (Lefty)
				lefty = grid.elem_find_neighbor(elem_num, 4);
				neighbor_num = lefty;}

			
			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n3);
				X2 = (*neighbor.n4);
				X3 = (*neighbor.n7);
				X4 = (*neighbor.n8);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);
					grid.natcoord(q, inter, elem_num);
					if (q.x >=  1.0) q.x =  0.99;
					if (q.x <= -1.0) q.x = -0.99;
					if (q.y >=  1.0) q.y =  0.99;
					if (q.y <= -1.0) q.y = -0.99;
					if (q.z >=  1.0) q.z =  0.99;
					if (q.z <= -1.0) q.z = -0.99;
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);

					return true;}}}}

	// Face 4
	if (face == 4){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Fronty's Toppy)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 5);}
			if (i == 2){		
				// Check Neighbor Face 2 (Toppy)
				toppy = grid.elem_find_neighbor(elem_num, 5);
				neighbor_num = toppy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Backy's Toppy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 5);}
			if (i == 4){		
				// Check Neighbor Face 4 (Backy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = backy;}
			if (i == 5){		
				// Check Neighbor Face 5 (Backy's Bottomy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 6);}
			if (i == 6){		
				// Check Neighbor Face 6 (Bottomy)
				bottomy = grid.elem_find_neighbor(elem_num, 6);
				neighbor_num = bottomy;}	
			if (i == 7){		
				// Check Neighbor Face 7 (Fronty's Bottomy)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 6);}	
			if (i == 8){		
				// Check Neighbor Face 8 (Fronty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = fronty;}
			

			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n1);
				X2 = (*neighbor.n3);
				X3 = (*neighbor.n5);
				X4 = (*neighbor.n7);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);
					grid.natcoord(q, inter, elem_num);
					if (q.x >=  1.0) q.x =  0.99;
					if (q.x <= -1.0) q.x = -0.99;
					if (q.y >=  1.0) q.y =  0.99;
					if (q.y <= -1.0) q.y = -0.99;
					if (q.z >=  1.0) q.z =  0.99;
					if (q.z <= -1.0) q.z = -0.99;
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);

					return true;}}}}

	// Face 5
	if (face == 5){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Backy's Lefty)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 4);}
			if (i == 2){		
				// Check Neighbor Face 2 (Backy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = backy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Backy's Righty)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 2);}
			if (i == 4){		
				// Check Neighbor Face 4 (Righty)
				righty = grid.elem_find_neighbor(elem_num, 2);
				neighbor_num = righty;}
			if (i == 5){		
				// Check Neighbor Face 5 (Fronty's Righty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 2);}
			if (i == 6){		
				// Check Neighbor Face 6 (Fronty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = fronty;}	
			if (i == 7){		
				// Check Neighbor Face 7 (Fronty's Lefty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 4);}	
			if (i == 8){		
				// Check Neighbor Face 8 (Lefty)
				lefty = grid.elem_find_neighbor(elem_num, 4);
				neighbor_num = lefty;}
			

			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n5);
				X2 = (*neighbor.n6);
				X3 = (*neighbor.n7);
				X4 = (*neighbor.n8);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);
					grid.natcoord(q, inter, elem_num);
					if (q.x >=  1.0) q.x =  0.99;
					if (q.x <= -1.0) q.x = -0.99;
					if (q.y >=  1.0) q.y =  0.99;
					if (q.y <= -1.0) q.y = -0.99;
					if (q.z >=  1.0) q.z =  0.99;
					if (q.z <= -1.0) q.z = -0.99;
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);

					return true;}}}}

	// Face 6
	if (face == 6){
		for (int i = 1; i < 9; i++){
			if (i == 1){		
				// Check Neighbor Face 1 (Backy's Lefty)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 4);}
			if (i == 2){		
				// Check Neighbor Face 2 (Backy)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = backy;}	
			if (i == 3){		
				// Check Neighbor Face 3 (Backy's Righty)
				backy = grid.elem_find_neighbor(elem_num, 3);
				neighbor_num = grid.elem_find_neighbor(backy, 2);}
			if (i == 4){		
				// Check Neighbor Face 4 (Righty)
				righty = grid.elem_find_neighbor(elem_num, 2);
				neighbor_num = righty;}
			if (i == 5){		
				// Check Neighbor Face 5 (Fronty's Righty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 2);}
			if (i == 6){		
				// Check Neighbor Face 6 (Fronty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = fronty;}	
			if (i == 7){		
				// Check Neighbor Face 7 (Fronty's Lefty)
				fronty = grid.elem_find_neighbor(elem_num, 1);
				neighbor_num = grid.elem_find_neighbor(fronty, 4);}	
			if (i == 8){		
				// Check Neighbor Face 8 (Lefty)
				lefty = grid.elem_find_neighbor(elem_num, 4);
				neighbor_num = lefty;}
			

			if (neighbor_num != -1){
				neighbor = grid.GetElement(neighbor_num);

				X1 = (*neighbor.n1);
				X2 = (*neighbor.n2);
				X3 = (*neighbor.n3);
				X4 = (*neighbor.n4);

				newton_find_intersect(new_lam, new_e1, new_e2, A, B, X1, X2, X3, X4);

				if (((new_lam > 0.0) && (new_lam <= 1.0)) && (fabs(new_e1) <= 1.0) && (fabs(new_e2) <= 1.0)){					
					lam = new_lam;
					e1 = new_e1;
					e2 = new_e2;
					inter.x = shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x;
					inter.y = shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y;
					inter.z = shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z;

					vec3d q(0,0,0);
					grid.natcoord(q, inter, elem_num);
					if (q.x >=  1.0) q.x =  0.99;
					if (q.x <= -1.0) q.x = -0.99;
					if (q.y >=  1.0) q.y =  0.99;
					if (q.y <= -1.0) q.y = -0.99;
					if (q.z >=  1.0) q.z =  0.99;
					if (q.z <= -1.0) q.z = -0.99;
					grid.nattoglobal(inter, q, elem_num);
					elem = grid.GetElement(neighbor_num);
					return true;}}}}

	return false;
	
}





///////////////////////////////////////////////////////////////////////
// newton_find_instersect
///////////////////////////////////////////////////////////////////////

void BC::newton_find_intersect(double &lam, double &e1, double &e2, vec3d &A, vec3d &B, Node &X1, Node &X2, Node &X3, Node &X4) 
{
    double err = 1;
	double tol = 1e-5;
	int iter = 0;
	int max_iter = 5;
	
	mat3d Hmat;
    vec3d R;
	vec3d E;
    vec3d dE;
    vec3d newE;
	
	E.x = lam;	// Initial lambda
	E.y = e1;  // Initial e1
	E.z = e2;  // Initial e2
	

	while ((err > tol) && (iter < max_iter)){
		lam = E.x;
		e1 = E.y;
		e2 = E.z;

		R.x = A.x + lam*(B.x - A.x) - ( shape_2D(1, e1, e2)*X1.rt.x + shape_2D(2, e1, e2)*X2.rt.x + shape_2D(3, e1, e2)*X3.rt.x + shape_2D(4, e1, e2)*X4.rt.x );
		R.y = A.y + lam*(B.y - A.y) - ( shape_2D(1, e1, e2)*X1.rt.y + shape_2D(2, e1, e2)*X2.rt.y + shape_2D(3, e1, e2)*X3.rt.y + shape_2D(4, e1, e2)*X4.rt.y );
		R.z = A.z + lam*(B.z - A.z) - ( shape_2D(1, e1, e2)*X1.rt.z + shape_2D(2, e1, e2)*X2.rt.z + shape_2D(3, e1, e2)*X3.rt.z + shape_2D(4, e1, e2)*X4.rt.z );

		Hmat[0][0] = B.x - A.x;
		Hmat[0][1] = - ( d1_shape_2D(1, 1, e1, e2)*X1.rt.x + d1_shape_2D(2, 1, e1, e2)*X2.rt.x + d1_shape_2D(3, 1, e1, e2)*X3.rt.x + d1_shape_2D(4, 1, e1, e2)*X4.rt.x );
		Hmat[0][2] = - ( d1_shape_2D(1, 2, e1, e2)*X1.rt.x + d1_shape_2D(2, 2, e1, e2)*X2.rt.x + d1_shape_2D(3, 2, e1, e2)*X3.rt.x + d1_shape_2D(4, 2, e1, e2)*X4.rt.x );

		Hmat[1][0] = B.y - A.y;
		Hmat[1][1] = - ( d1_shape_2D(1, 1, e1, e2)*X1.rt.y + d1_shape_2D(2, 1, e1, e2)*X2.rt.y + d1_shape_2D(3, 1, e1, e2)*X3.rt.y + d1_shape_2D(4, 1, e1, e2)*X4.rt.y );
		Hmat[1][2] = - ( d1_shape_2D(1, 2, e1, e2)*X1.rt.y + d1_shape_2D(2, 2, e1, e2)*X2.rt.y + d1_shape_2D(3, 2, e1, e2)*X3.rt.y + d1_shape_2D(4, 2, e1, e2)*X4.rt.y );

		Hmat[2][0] = B.z - A.z;
		Hmat[2][1] = - ( d1_shape_2D(1, 1, e1, e2)*X1.rt.z + d1_shape_2D(2, 1, e1, e2)*X2.rt.z + d1_shape_2D(3, 1, e1, e2)*X3.rt.z + d1_shape_2D(4, 1, e1, e2)*X4.rt.z );
		Hmat[2][2] = - ( d1_shape_2D(1, 2, e1, e2)*X1.rt.z + d1_shape_2D(2, 2, e1, e2)*X2.rt.z + d1_shape_2D(3, 2, e1, e2)*X3.rt.z + d1_shape_2D(4, 2, e1, e2)*X4.rt.z );

		dE = Hmat.inverse()*R;
		newE = E - dE;
        
		err = dE.norm();
		E = newE;
		++iter;}

	return;
}



///////////////////////////////////////////////////////////////////////
// shape_2D
///////////////////////////////////////////////////////////////////////

double BC::shape_2D(int node, double e1, double e2)
{
	double N = 0.;

	if (node == 1){
		N = (1./4.)*(1 - e1)*(1 - e2);
		return N;}

	if (node == 2){
		N = (1./4.)*(1 + e1)*(1 - e2);
		return N;}

	if (node == 3){
		N = (1./4.)*(1 - e1)*(1 + e2);
		return N;}

	if (node == 4){
		N = (1./4.)*(1 + e1)*(1 + e2);
		return N;}

	return N;
}



///////////////////////////////////////////////////////////////////////
// d1_shape_2D
///////////////////////////////////////////////////////////////////////

double BC::d1_shape_2D(int node, int d, double e1, double e2)
{
	double dN = 0.;

	if ((node == 1) && (d == 1)){
		dN = -(1./4.)*(1 - e2);
		return dN;}
	if ((node == 1) && (d == 2)){
		dN = -(1./4.)*(1. - e1);
		return dN;}
	
	if ((node == 2) && (d == 1)){
		dN = (1./4.)*(1. - e2);
		return dN;}
	if ((node == 2) && (d == 2)){
		dN = -(1./4.)*(1. + e1);
		return dN;}

	if ((node == 3) && (d == 1)){
		dN = -(1./4.)*(1 + e2);
		return dN;}
	if ((node == 3) && (d == 2)){
		dN = (1./4.)*(1 - e1);
		return dN;}

	if ((node == 4) && (d == 1)){
		dN = (1./4.)*(1 + e2);
		return dN;}
	if ((node == 4) && (d == 2)){
		dN = (1./4.)*(1 + e1);
		return dN;}

	return dN;
}



///////////////////////////////////////////////////////////////////////
// intersceptface (old method, not used currently)
///////////////////////////////////////////////////////////////////////

vec3d BC::intersceptface(int face, double &xix_0, double &xiy_0, double &xiz_0, double &xix_1, double &xiy_1, double &xiz_1, Segment &seg, int k)
{
    vec3d i_point;
    double xix_i, xiy_i, xiz_i = {0.};
    
    vec3d p;
    vec3d n;
    vec3d v;
    vec3d v2;
    vec3d lp;
    double u;
    
    if (k == 1){
        v.x = xix_1 - xix_0;
        v.y = xiy_1 - xiy_0;
        v.z = xiz_1 - xiy_0;
        lp.x = xix_0;
        lp.y = xiy_0;
        lp.z = xiz_0;}
    else if (k == 0){
	    v.x = xix_0 - xix_1;
        v.y = xiy_0 - xiy_1;
        v.z = xiz_0 - xiy_1;
        lp.x = xix_1;
        lp.y = xiy_1;
        lp.z = xiz_1;}
    
    v = v/v.norm();
        
    if (face == 1){
        n.x = 0;
        n.y = -1;
        n.z = 0;
        p.x = 0;
        p.y = -1;
        p.z = 0;}
        
    if (face == 2){
        n.x = 1;
        n.y = 0;
        n.z = 0;
        p.x = 1;
        p.y = 0;
        p.z = 0;}
        
    if (face == 3){
        n.x = 0;
        n.y = 1;
        n.z = 0;
        p.x = 0;
        p.y = 1;
        p.z = 0;}
                
    if (face == 4){
        n.x = -1;
        n.y = 0;
        n.z = 0;
        p.x = -1;
        p.y = 0;
        p.z = 0;}
                    
    if (face == 5){
        n.x = 0;
        n.y = 0;
        n.z = 1;
        p.x = 0;
        p.y = 0;
        p.z = 1;}        
    
    if (face == 6){
        n.x = 0;
        n.y = 0;
        n.z = -1;
        p.x = 0;
        p.y = 0;
        p.z = -1;}
      
    
    v2 = p - lp;
    
    u = (n*v2) / (n*v);
    
    
	if (u == 0){
	    u = 0.01;
	    seg.tip(0).bactive = false;
	    seg.tip(1).bactive = false;
	}
	
	xix_i = lp.x + u*v.x;
	xiy_i = lp.y + u*v.y;
	xiz_i = lp.z + u*v.z;
	                    
    if (xix_i > 1.)
        xix_i = 1.;
    else if (xix_i < -1.)
        xix_i = -1.;
    
    if (xiy_i > 1.)
        xiy_i = 1.;
    else if (xiy_i < -1.)
        xiy_i = -1.; 
    
    if (xiz_i > 1.)
        xiz_i = 1.;
    else if (xiz_i < -1.)
        xiz_i = -1.;
    
    i_point.x = xix_i;
    i_point.y = xiy_i;
    i_point.z = xiz_i;
    
    return i_point;
}


///////////////////////////////////////////////////////////////////////
// collfibwallBC
///////////////////////////////////////////////////////////////////////

void BC::collfibwallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k)
{
 //   double xpt_0, ypt_0, zpt_0 = {0.};
 //   double xpt_1, ypt_1, zpt_1 = {0.};      
	//double xpt_i = 0.; double ypt_i = 0.; double zpt_i = 0.;

 //   xpt_0 = seg.x[0];
 //   ypt_0 = seg.y[0];
 //   zpt_0 = seg.z[0];
 //   xpt_1 = seg.x[1];
 //   ypt_1 = seg.y[1];
 //   zpt_1 = seg.z[1];
 //
	//xpt_i = i_point.x;
	//ypt_i = i_point.y;
	//zpt_i = i_point.z;

 //   Segment seg2;
 //           
 //   if (k == 1){
 //       seg.x[1] = xpt_i;
 //       seg.y[1] = ypt_i;
 //       seg.z[1] = zpt_i;
 //       seg.findlength();
 //       seg.tip[1] = 0;
 //       seg2.tip[1] = 1;
	//	seg2.tip_elem[0] = seg.tip_elem[1];}
 //   else if (k == 0){
 //       seg.x[0] = xpt_i;
 //       seg.y[0] = ypt_i;
 //       seg.z[0] = zpt_i;
 //       seg.findlength();
 //       seg.tip[0] = 0;
 //       seg2.tip[0] = -1;
	//	seg2.tip_elem[1] = seg.tip_elem[0];}

 //   seg.BCdead = 1;
 //      
 //   if (k == 1){
 //       seg2.x[0] = xpt_i;
 //       seg2.y[0] = ypt_i;
 //       seg2.z[0] = zpt_i;
 //       seg2.x[1] = xpt_i + (xpt_1 - xpt_i);
 //       seg2.y[1] = ypt_i + (ypt_1 - ypt_i);
 //       seg2.z[1] = zpt_i + (zpt_1 - zpt_i);}
 //   else if (k == 0){
 //       seg2.x[1] = xpt_i;
 //       seg2.y[1] = ypt_i;
 //       seg2.z[1] = zpt_i;
 //       seg2.x[0] = xpt_i + (xpt_0 - xpt_i);
 //       seg2.y[0] = ypt_i + (ypt_0 - ypt_i);
 //       seg2.z[0] = zpt_i + (zpt_0 - zpt_i);}
 //   
 //   
 //   switch (face)
 //   {
 //   case 1:
 //       if (k == 1)
 //           seg2.y[1] = ypt_i - (ypt_1 - ypt_i);
 //       else if (k == 0)
 //           seg2.y[0] = ypt_i - (ypt_0 - ypt_i);
 //       break;
 //   
 //   case 2:
 //       if (k == 1)
 //           seg2.x[1] = xpt_i - (xpt_1 - xpt_i);
 //       else if (k == 0)
 //           seg2.x[0] = xpt_i - (xpt_0 - xpt_i);
 //       break;
 //           
 //   case 3:
 //       if (k == 1)
 //           seg2.y[1] = ypt_i - (ypt_1 - ypt_i);
 //       else if (k == 0)
 //           seg2.y[0] = ypt_i - (ypt_0 - ypt_i);
 //       break;       
 //   
 //   case 4:
 //       if (k == 1)
 //           seg2.x[1] = xpt_i - (xpt_1 - xpt_i);
 //       else if (k == 0)
 //           seg2.x[0] = xpt_i - (xpt_0 - xpt_i);
 //       break;
 //   
 //   case 5:
 //       if (k == 1)
 //           seg2.z[1] = zpt_i - (zpt_1 - zpt_i);
 //       else if (k == 0)
 //           seg2.z[0] = zpt_i - (zpt_0 - zpt_i);
 //       break;     
 //   
 //   case 6:
 //       if (k == 1)
 //           seg2.z[1] = zpt_i - (zpt_1 - zpt_i);
 //       else if (k == 0)
 //           seg2.z[0] = zpt_i - (zpt_0 - zpt_i);
 //       break;
 //   
	//case -1:
	//	seg2.tip[0] = 0;
	//	seg2.tip[1] = 0;
	//	seg2.mark_of_death = true;
	//	seg2.death_label = 6;
	//	return;
	//	break;
	//}
 //              
 //   seg2.BCdead = 1;
 //   seg2.m_nseed = seg.m_nseed;
 //   seg2.TofBirth = data.t;
 //   data.num_vessel++; 
 //   seg2.sprout = seg.sprout;
 //   
 //   seg2.findlength();
 //   seg2.findphi();
 //   seg2.bdyf_id[k] = seg.bdyf_id[k]; 

	//elem_num = grid.findelem(seg2.x[k], seg2.y[k], seg2.z[k]);
 //       
 //   Culture cult;
 //   
 //   double coll_phi1, coll_phi2 = {0.};
 //   double W1 = 10*(1/grid.den_scale);
 //   double W2 = 100;
 //   
 //   if (k == 1){
 //       coll_phi1 = cult.findCollAngle(seg2.x[0], seg2.y[0], seg2.z[0], grid, data, 1);
 //       coll_phi2 = cult.findCollAngle(seg2.x[0], seg2.y[0], seg2.z[0], grid, data, 2);
 //       //seg2.phi1 = W1*coll_phi1 + W2*seg2.phi1;
 //       //seg2.phi2 = W1*coll_phi2 + W2*seg2.phi2;
 //       seg2.phi1 = coll_phi1;
 //       seg2.phi2 = coll_phi2;
 //       seg2.x[1] = seg2.x[0] + seg2.length*cos(seg2.phi2)*cos(seg2.phi1);
 //       seg2.y[1] = seg2.y[0] + seg2.length*cos(seg2.phi2)*sin(seg2.phi1);
 //       seg2.z[1] = seg2.z[0] + seg2.length*sin(seg2.phi2);}
 //   else if (k == 0){
 //       coll_phi1 = cult.findCollAngle(seg2.x[1], seg2.y[1], seg2.z[1], grid, data, 1);
 //       coll_phi2 = cult.findCollAngle(seg2.x[1], seg2.y[1], seg2.z[1], grid, data, 2);
 //       /*seg2.phi1 = W1*coll_phi1 + W2*seg2.phi1;
 //       seg2.phi2 = W1*coll_phi2 + W2*seg2.phi2;*/
 //       seg2.phi1 = coll_phi1;
 //       seg2.phi2 = coll_phi2;
 //       seg2.x[0] = seg2.x[1] + seg2.length*cos(seg2.phi2)*cos(seg2.phi1);
 //       seg2.y[0] = seg2.y[1] + seg2.length*cos(seg2.phi2)*sin(seg2.phi1);
 //       seg2.z[0] = seg2.z[1] + seg2.length*sin(seg2.phi2);}
                  
    return;
}



///////////////////////////////////////////////////////////////////////
// inplanewallBC
///////////////////////////////////////////////////////////////////////

Segment BC::inplanewallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k)
{
	Culture& cult = m_angio.GetCulture();

    BC_bouncy = true;
	
	double xpt_0, ypt_0, zpt_0 = {0.};
    double xpt_1, ypt_1, zpt_1 = {0.};      
	double xpt_i = 0.; double ypt_i = 0.; double zpt_i = 0.;

	double old_length = 0.;
	double remain_length = 0.;

    xpt_0 = seg.tip(0).rt.x;
    ypt_0 = seg.tip(0).rt.y;
    zpt_0 = seg.tip(0).rt.z;
    xpt_1 = seg.tip(1).rt.x;
    ypt_1 = seg.tip(1).rt.y;
    zpt_1 = seg.tip(1).rt.z;
 
	xpt_i = i_point.x;
	ypt_i = i_point.y;
	zpt_i = i_point.z;

	double eps = 1e-3;

	if (face == 1){
		ypt_i = (1.0 + eps)*i_point.y;}
	if (face == 2){
		xpt_i = (1.0 - eps)*i_point.x;}
	if (face == 3){
		ypt_i = (1.0 - eps)*i_point.y;}
	if (face == 4){
		xpt_i = (1.0 + eps)*i_point.x;}
	if (face == 5){
		zpt_i = (1.0 - eps)*i_point.z;}
	if (face == 6){
		zpt_i = (1.0 + eps)*i_point.z;}
		
	old_length = seg.length();
	
	Segment seg2;
    seg2.m_nseed = seg.m_nseed;

    if (k == 1){
        seg.tip(1).rt.x = xpt_i;
        seg.tip(1).rt.y = ypt_i;
        seg.tip(1).rt.z = zpt_i;
        seg.Update();
        seg.tip(1).bactive = false;
        seg2.tip(1).bactive = true;
		seg2.tip(0).pt = seg.tip(1).pt;
	}
    else if (k == 0){
        seg.tip(0).rt.x = xpt_i;
        seg.tip(0).rt.y = ypt_i;
        seg.tip(0).rt.z = zpt_i;
        seg.Update();
        seg.tip(0).bactive = false;
        seg2.tip(0).bactive = true;
		seg2.tip(1).pt = seg.tip(0).pt;
	}

    seg.SetFlagOn(Segment::BC_DEAD);
       
    remain_length = fabs(old_length) - fabs(seg.length());

	if (k == 0)
		remain_length = -1.0*remain_length;
	
	vec3d face_normal;
	vec3d s1; vec3d s2; vec3d s3; vec3d s4;
	vec3d c1; vec3d c2; vec3d c3; vec3d c4;

	Grid& grid = m_angio.GetGrid();
	Elem elem;
	elem = grid.GetElement(elem_num);
	
	switch (face)
    {
    case 1:
		s1.x = (*elem.n5).rt.x - (*elem.n1).rt.x; s1.y = (*elem.n5).rt.y - (*elem.n1).rt.y; s1.z = (*elem.n5).rt.z - (*elem.n1).rt.z;     
		s2.x = (*elem.n6).rt.x - (*elem.n5).rt.x; s2.y = (*elem.n6).rt.y - (*elem.n5).rt.y; s2.z = (*elem.n6).rt.z - (*elem.n5).rt.z;  
		s3.x = (*elem.n2).rt.x - (*elem.n6).rt.x; s3.y = (*elem.n2).rt.y - (*elem.n6).rt.y; s3.z = (*elem.n2).rt.z - (*elem.n6).rt.z;  
		s4.x = (*elem.n1).rt.x - (*elem.n2).rt.x; s4.y = (*elem.n1).rt.y - (*elem.n2).rt.y; s4.z = (*elem.n1).rt.z - (*elem.n2).rt.z;  
        break;
    
    case 2:
        s1.x = (*elem.n8).rt.x - (*elem.n6).rt.x; s1.y = (*elem.n8).rt.y - (*elem.n6).rt.y; s1.z = (*elem.n8).rt.z - (*elem.n6).rt.z;     
		s2.x = (*elem.n4).rt.x - (*elem.n8).rt.x; s2.y = (*elem.n4).rt.y - (*elem.n8).rt.y; s2.z = (*elem.n4).rt.z - (*elem.n8).rt.z;  
		s3.x = (*elem.n2).rt.x - (*elem.n4).rt.x; s3.y = (*elem.n2).rt.y - (*elem.n4).rt.y; s3.z = (*elem.n2).rt.z - (*elem.n4).rt.z;  
		s4.x = (*elem.n6).rt.x - (*elem.n2).rt.x; s4.y = (*elem.n6).rt.y - (*elem.n2).rt.y; s4.z = (*elem.n6).rt.z - (*elem.n2).rt.z;  
        break;
            
    case 3:
        s1.x = (*elem.n8).rt.x - (*elem.n4).rt.x; s1.y = (*elem.n8).rt.y - (*elem.n4).rt.y; s1.z = (*elem.n8).rt.z - (*elem.n4).rt.z;     
		s2.x = (*elem.n7).rt.x - (*elem.n8).rt.x; s2.y = (*elem.n7).rt.y - (*elem.n8).rt.y; s2.z = (*elem.n7).rt.z - (*elem.n8).rt.z;  
		s3.x = (*elem.n3).rt.x - (*elem.n7).rt.x; s3.y = (*elem.n3).rt.y - (*elem.n7).rt.y; s3.z = (*elem.n3).rt.z - (*elem.n7).rt.z;  
		s4.x = (*elem.n4).rt.x - (*elem.n3).rt.x; s4.y = (*elem.n4).rt.y - (*elem.n3).rt.y; s4.z = (*elem.n4).rt.z - (*elem.n3).rt.z;  
        break;       
    
    case 4:
        s1.x = (*elem.n7).rt.x - (*elem.n3).rt.x; s1.y = (*elem.n7).rt.y - (*elem.n3).rt.y; s1.z = (*elem.n7).rt.z - (*elem.n3).rt.z;     
		s2.x = (*elem.n5).rt.x - (*elem.n7).rt.x; s2.y = (*elem.n5).rt.y - (*elem.n7).rt.y; s2.z = (*elem.n5).rt.z - (*elem.n7).rt.z;  
		s3.x = (*elem.n1).rt.x - (*elem.n5).rt.x; s3.y = (*elem.n1).rt.y - (*elem.n5).rt.y; s3.z = (*elem.n1).rt.z - (*elem.n5).rt.z;  
		s4.x = (*elem.n3).rt.x - (*elem.n1).rt.x; s4.y = (*elem.n3).rt.y - (*elem.n1).rt.y; s4.z = (*elem.n3).rt.z - (*elem.n1).rt.z;
        break;
    
    case 5:
        s1.x = (*elem.n7).rt.x - (*elem.n5).rt.x; s1.y = (*elem.n7).rt.y - (*elem.n5).rt.y; s1.z = (*elem.n7).rt.z - (*elem.n5).rt.z;     
		s2.x = (*elem.n8).rt.x - (*elem.n7).rt.x; s2.y = (*elem.n8).rt.y - (*elem.n7).rt.y; s2.z = (*elem.n8).rt.z - (*elem.n7).rt.z;  
		s3.x = (*elem.n6).rt.x - (*elem.n8).rt.x; s3.y = (*elem.n6).rt.y - (*elem.n8).rt.y; s3.z = (*elem.n6).rt.z - (*elem.n8).rt.z;  
		s4.x = (*elem.n5).rt.x - (*elem.n6).rt.x; s4.y = (*elem.n5).rt.y - (*elem.n6).rt.y; s4.z = (*elem.n5).rt.z - (*elem.n6).rt.z;
        break;     
    
    case 6:
        s1.x = (*elem.n1).rt.x - (*elem.n3).rt.x; s1.y = (*elem.n1).rt.y - (*elem.n3).rt.y; s1.z = (*elem.n1).rt.z - (*elem.n3).rt.z;     
		s2.x = (*elem.n2).rt.x - (*elem.n1).rt.x; s2.y = (*elem.n2).rt.y - (*elem.n1).rt.y; s2.z = (*elem.n2).rt.z - (*elem.n1).rt.z;  
		s3.x = (*elem.n4).rt.x - (*elem.n2).rt.x; s3.y = (*elem.n4).rt.y - (*elem.n2).rt.y; s3.z = (*elem.n4).rt.z - (*elem.n2).rt.z;  
		s4.x = (*elem.n3).rt.x - (*elem.n4).rt.x; s4.y = (*elem.n3).rt.y - (*elem.n4).rt.y; s4.z = (*elem.n3).rt.z - (*elem.n4).rt.z;
        break;
    
	case -1:
		seg2.tip(0).bactive = false;
		seg2.tip(1).bactive = false;
		seg2.mark_of_death = true;
		seg2.death_label = 1;
		assert(false);
		return seg2;
		break;
	}
    
	s1 = s1/s1.norm();
	s2 = s2/s2.norm();
	s3 = s3/s3.norm();
	s4 = s4/s4.norm();
	
	c1 = s1^s4;
	c2 = s2^s1;
	c3 = s3^s2;
	c4 = s4^s3;
	
	face_normal = (c1 + c2 + c3 + c4)*0.25;

	vec3d seg_unit_perp; vec3d seg_unit_inplane;

	seg_unit_perp = face_normal*(seg.uvect()*face_normal);
	seg_unit_inplane = seg.uvect() - seg_unit_perp;
	
	if (k == 1){
        seg2.tip(0).rt.x = xpt_i;
        seg2.tip(0).rt.y = ypt_i;
        seg2.tip(0).rt.z = zpt_i;
		seg2.tip(1).rt.x = xpt_i + remain_length*seg_unit_inplane.x;
        seg2.tip(1).rt.y = ypt_i + remain_length*seg_unit_inplane.y;
        seg2.tip(1).rt.z = zpt_i + remain_length*seg_unit_inplane.z;}
    else if (k == 0){
        seg2.tip(1).rt.x = xpt_i;
        seg2.tip(1).rt.y = ypt_i;
        seg2.tip(1).rt.z = zpt_i;
        seg2.tip(0).rt.x = xpt_i + remain_length*seg_unit_inplane.x;
        seg2.tip(0).rt.y = ypt_i + remain_length*seg_unit_inplane.y;
        seg2.tip(0).rt.z = zpt_i + remain_length*seg_unit_inplane.z;}

//	seg2.m_length = remain_length;
	seg2.SetFlagOn(Segment::BC_DEAD);
    /*data.num_vessel++; 
	seg2.line_num = seg.line_num;*/

	// TODO: I don't think I want to copy all the flags. Maybe I don't need to copy anything.
	if (seg.GetFlag(Segment::INIT_SPROUT) == false)
		seg2.SetFlagOn(seg.GetFlags());
    
	seg2.Update();
	//seg2.findphi();
    seg2.tip(k).bdyf_id = seg.tip(k).bdyf_id; 

	return seg2;
}


///////////////////////////////////////////////////////////////////////
// inplanewallBC
///////////////////////////////////////////////////////////////////////

Segment BC::symplaneperiodicwallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k)
{
	Culture& cult = m_angio.GetCulture();

    BC_bouncy = true;
	
	double xpt_0, ypt_0, zpt_0 = {0.};
    double xpt_1, ypt_1, zpt_1 = {0.};      
	double xpt_i = 0.; double ypt_i = 0.; double zpt_i = 0.;

	double old_length = 0.;
	double remain_length = 0.;

    xpt_0 = seg.tip(0).rt.x;
    ypt_0 = seg.tip(0).rt.y;
    zpt_0 = seg.tip(0).rt.z;
    xpt_1 = seg.tip(1).rt.x;
    ypt_1 = seg.tip(1).rt.y;
    zpt_1 = seg.tip(1).rt.z;
 
	xpt_i = i_point.x;
	ypt_i = i_point.y;
	zpt_i = i_point.z;

	double eps = 1e-3;

	if (face == 1){
		ypt_i = (1.0 + eps)*i_point.y;}
	if (face == 2){
		xpt_i = (1.0 - eps)*i_point.x;}
	if (face == 3){
		ypt_i = (1.0 - eps)*i_point.y;}
	if (face == 4){
		xpt_i = (1.0 + eps)*i_point.x;}
	if (face == 5){
		zpt_i = (1.0 - eps)*i_point.z;}
	if (face == 6){
		zpt_i = (1.0 + eps)*i_point.z;}
		
	old_length = seg.length();
	
	Segment seg2;
    seg2.m_nseed = seg.m_nseed;

    if (k == 1){
        seg.tip(1).rt.x = xpt_i;
        seg.tip(1).rt.y = ypt_i;
        seg.tip(1).rt.z = zpt_i;
        seg.Update();
        seg.tip(1).bactive = false;
		seg2.tip(0).pt = seg.tip(1).pt;
	}
    else if (k == 0){
        seg.tip(0).rt.x = xpt_i;
        seg.tip(0).rt.y = ypt_i;
        seg.tip(0).rt.z = zpt_i;
        seg.Update();
        seg2.tip(0).bactive = true;
		seg2.tip(1).pt = seg.tip(0).pt;
	}

    seg.SetFlagOn(Segment::BC_DEAD);
       
    remain_length = fabs(old_length) - seg.length();

	vec3d new_uvect;

	new_uvect = seg.uvect();
	
	double new_xpt; double new_ypt; double new_zpt;
	vec3d rand_disp;
	double disp_dist_min = 50;
	double disp_dist_max = 100;
	double disp_dist;

	rand_disp = vrand();
	disp_dist = disp_dist_min + frand()*(disp_dist_max - disp_dist_min);
	
	if ((face == 1) || (face == 3)){
		new_uvect.y = -1.0*new_uvect.y;
		rand_disp.y = 0.;}

	if ((face == 2) || (face == 4)){
		new_uvect.x = -1.0*new_uvect.x;
		rand_disp.x = 0.;}

	if ((face == 5) || (face == 6)){
		new_uvect.z = -1.0*new_uvect.z;
		rand_disp.z = 0.;}

	rand_disp.x = disp_dist*rand_disp.x; 
	rand_disp.y = disp_dist*rand_disp.y; 
	rand_disp.z = disp_dist*rand_disp.z;
	
	new_xpt = xpt_i + rand_disp.x; 
	new_ypt = ypt_i + rand_disp.y; 
	new_zpt = zpt_i + rand_disp.z;

    seg2.tip(0).rt.x = new_xpt;
    seg2.tip(0).rt.y = new_ypt;
    seg2.tip(0).rt.z = new_zpt;
	seg2.tip(1).rt.x = new_xpt + remain_length*new_uvect.x;
    seg2.tip(1).rt.y = new_ypt + remain_length*new_uvect.y;
    seg2.tip(1).rt.z = new_zpt + remain_length*new_uvect.z;

	seg2.tip(1).bactive = true;
//	seg2.m_length = remain_length;
	seg2.SetFlagOn(Segment::BC_DEAD);

	// TODO: I don't think I want to copy all the flags. Maybe I don't need to copy anything.
	if (seg.GetFlag(Segment::INIT_SPROUT) == false)
		seg2.SetFlagOn(seg.GetFlags());
    
	seg2.Update();
    seg2.tip(k).bdyf_id = seg.tip(k).bdyf_id; 

	return seg2;
}


///////////////////////////////////////////////////////////////////////
// PeriodicBC
///////////////////////////////////////////////////////////////////////

//table of faces opposite to 0,1,..6
double oppface[6];
// TODO: vessel lenghts are always positive. Fix the logic here.
Segment BC::PeriodicBC(Segment &seg)
{
	assert(false);
/*	Grid& grid = m_angio.GetGrid();
	Culture& cult = m_angio.GetCulture();

	oppface[0] = grid.yrange[1];
	oppface[1] = grid.xrange[0];
	oppface[2] = grid.yrange[0];
	oppface[3] = grid.xrange[1];
	oppface[4] = grid.zrange[0];
	oppface[5] = grid.zrange[1];
	double unit_vec[3] = {0};
	double length = 0.0;
	double rem_length = 0.0;
	int n = 0;
	double intersectpt[3] = {0};
	
	if (seg.tip(1).bactive)
	{
		unit_vec[0] = (seg.tip(1).rt.x-seg.tip(0).rt.x);
		unit_vec[1] = (seg.tip(1).rt.y-seg.tip(0).rt.y);
		unit_vec[2] = (seg.tip(1).rt.z-seg.tip(0).rt.z);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		
		for (n=0;n<6;++n)
		{
			if (grid.intersectPlane(seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.m_nseed << " Crossed Face " << n << endl;
				seg.tip(1).rt.x = intersectpt[0];
				seg.tip(1).rt.y = intersectpt[1];
				seg.tip(1).rt.z = intersectpt[2];
				seg.Update();
				if (seg.length() < 0)
				{
					seg.m_length = -(seg.tip(1).rt - seg.tip(0).rt).norm();
				}
				else
				{
					seg.m_length = (seg.tip(1).rt - seg.tip(0).rt).norm();
				}

				seg.tip(1).bactive = false;
				seg.SetFlagOn(Segment::BC_DEAD);
				
				if (n > 3)
				    cult.m_num_zdead += 1;
				return seg;
				
				cult.AddSegment(seg);
				Segment seg2;
				cult.m_num_vessel += 1;
				
				switch (n)
				{
				case 0:
					seg2.tip(0).rt.x = seg.tip(1).rt.x;
					seg2.tip(0).rt.y = oppface[n];
					seg2.tip(0).rt.z = seg.tip(1).rt.z;
					break;
					
				case 1:
					seg2.tip(0).rt.x = oppface[n];
					seg2.tip(0).rt.y = seg.tip(1).rt.y;
					seg2.tip(0).rt.z = seg.tip(1).rt.z;
					break;
				case 2:
					seg2.tip(0).rt.x = seg.tip(1).rt.x;
					seg2.tip(0).rt.y = oppface[n];
					seg2.tip(0).rt.z = seg.tip(1).rt.z;
					break;
				case 3:
					seg2.tip(0).rt.x = oppface[n];
					seg2.tip(0).rt.y = seg.tip(1).rt.y;
					seg2.tip(0).rt.z = seg.tip(1).rt.z;
					break;
				case 4:
					seg2.tip(0).rt.x = seg.tip(1).rt.x;
					seg2.tip(0).rt.y = seg.tip(1).rt.y;
					seg2.tip(0).rt.z = oppface[n];
					break;
				case 5:
					seg2.tip(0).rt.x = seg.tip(1).rt.x;
					seg2.tip(0).rt.y = seg.tip(1).rt.y;
					seg2.tip(0).rt.z = oppface[n];
					break;
				}

				if (seg.length() > 0) 
				{
					rem_length = length - (seg.tip(1).rt - seg.tip(0).rt).norm();
					seg2.tip(1).rt.x = seg2.tip(0).rt.x + rem_length*unit_vec[0];
					seg2.tip(1).rt.y = seg2.tip(0).rt.y + rem_length*unit_vec[1];
					seg2.tip(1).rt.z = seg2.tip(0).rt.z + rem_length*unit_vec[2];
				}
				else 
				{
					rem_length = -length + (seg.tip(1).rt- seg.tip(0).rt).norm();
					seg2.tip(1).rt.x = seg2.tip(0).rt.x - rem_length*unit_vec[0];
					seg2.tip(1).rt.y = seg2.tip(0).rt.y - rem_length*unit_vec[1];
					seg2.tip(1).rt.z = seg2.tip(0).rt.z - rem_length*unit_vec[2];
				}
				
				seg2.tip(1).bactive = true;
				seg2.tip(0).bactive = false;
				seg2.m_nseed = seg.m_nseed;
				seg2.m_nvessel = cult.m_num_vessel;
				seg2.SetFlagOn(seg.GetFlags());
//				seg2.m_length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2);
				else
					return seg2;
			}
			
		}
	}
	
	else
	{
		unit_vec[0] = (seg.tip(0).rt.x-seg.tip(1).rt.x);
		unit_vec[1] = (seg.tip(0).rt.y-seg.tip(1).rt.y);
		unit_vec[2] = (seg.tip(0).rt.z-seg.tip(1).rt.z);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		for (n=0;n<6;++n)
		{
			if (grid.intersectPlane(seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.m_nseed << " Crossed Face " << n << endl;
				seg.tip(0).rt.x = intersectpt[0];
				seg.tip(0).rt.y = intersectpt[1];
				seg.tip(0).rt.z = intersectpt[2];
				if (seg.length() < 0)
				{
//					seg.m_length = -(seg.tip(1).rt - seg.tip(0).rt).norm();
				}
				else
				{
//					seg.m_length = (seg.tip(1).rt - seg.tip(0).rt).norm();
				}
				
				seg.tip(0).bactive = false;
                
				if (n > 3)
				    cult.m_num_zdead += 1;
				return seg;
				
				cult.AddSegment(seg);
				Segment seg2;
				cult.m_num_vessel += 1;
				
				switch (n)
				{
				case 0:
					seg2.tip(1).rt.x = seg.tip(0).rt.x;
					seg2.tip(1).rt.y = oppface[n];
					seg2.tip(1).rt.z = seg.tip(0).rt.z;
					break;
					
				case 1:
					seg2.tip(1).rt.x = oppface[n];
					seg2.tip(1).rt.y = seg.tip(0).rt.y;
					seg2.tip(1).rt.z = seg.tip(0).rt.z;
					break;
				case 2:
					seg2.tip(1).rt.x = seg.tip(0).rt.x;
					seg2.tip(1).rt.y = oppface[n];
					seg2.tip(1).rt.z = seg.tip(0).rt.z;
					break;
				case 3:
					seg2.tip(1).rt.x = oppface[n];
					seg2.tip(1).rt.y = seg.tip(0).rt.y;
					seg2.tip(1).rt.z = seg.tip(0).rt.z;
					break;
				case 4:
					seg2.tip(1).rt.x = seg.tip(0).rt.x;
					seg2.tip(1).rt.y = seg.tip(0).rt.y;
					seg2.tip(1).rt.z = oppface[n];
					break;
				case 5:
					seg2.tip(1).rt.x = seg.tip(0).rt.x;
					seg2.tip(1).rt.y = seg.tip(0).rt.y;
					seg2.tip(1).rt.z = oppface[n];
					break;
				}

				if (seg.length() < 0) 
				{
					rem_length = -length + (seg.tip(1).rt - seg.tip(0).rt).norm();
					seg2.tip(0).rt.x = seg2.tip(1).rt.x - rem_length*unit_vec[0];
					seg2.tip(0).rt.y = seg2.tip(1).rt.y - rem_length*unit_vec[1];
					seg2.tip(0).rt.z = seg2.tip(1).rt.z - rem_length*unit_vec[2];
					
				}
				else 
				{
					rem_length = length - (seg.tip(1).rt - seg.tip(0).rt).norm();
					seg2.tip(0).rt.x = seg2.tip(1).rt.x + rem_length*unit_vec[0];
					seg2.tip(0).rt.y = seg2.tip(1).rt.y + rem_length*unit_vec[1];
					seg2.tip(0).rt.z = seg2.tip(1).rt.z + rem_length*unit_vec[2];
				}
				
				seg2.tip(1).bactive = false;
				seg2.tip(0).bactive = true;
				seg2.m_nseed = seg.m_nseed;
				seg2.m_nvessel = cult.m_num_vessel;
				seg2.SetFlagOn(seg.GetFlags());
//				seg2.m_length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2);
				else
					return seg2;
			}
		}
	}
*/
	return seg;
}
