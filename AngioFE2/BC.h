#pragma once
#include <list>
using namespace std;
#include <FECore/vec3d.h>
#include "Grid.h"

//-----------------------------------------------------------------------------
class FEAngio;
class Elem;
class Segment;
class Node;

//-----------------------------------------------------------------------------
// The BC class is used by the CULTURE class to handle boundary 
// conditions within the model.
class BC
{
public:
	enum {
		STOP,		// the segment stops at the boundary
		BOUNCY		// bouncy wall
	};

public:
	BC(FEAngio& angio);
	~BC();
	
	// checks if a new segment has cross the boundary
	void CheckBC(Segment &seg);

private:
	// enforces a boundary condition
	void EnforceBC(Segment &seg, FACE_INTERSECTION& ic);

	// apply BC where vessel stops growing after hitting boundary
	void BCStop(Segment &seg, FACE_INTERSECTION& ic);

	// apply BC where the vessel bounces off the boundary
	void BCBouncy(Segment &seg, FACE_INTERSECTION& ic);

public:
	void collfibwallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k);
	vec3d intersceptface(int face, double &xix_0, double &xiy_0, double &xiz_0, double &xix_1, double &xiy_1, double &xiz_1, Segment &seg, int k);
	vec3d find_intersect(Elem &elem, int &face, Segment &seg);
	bool search_neighbors_4_intersect(Elem &elem, int face, double &lam, double &e1, double &e2, vec3d &A, vec3d &B, vec3d &inter); 
	void newton_find_intersect(double &lam, double &e1, double &e2, vec3d &A, vec3d &B, Node &X1, Node &X2, Node &X3, Node &X4);
	double shape_2D(int node, double e1, double e2);
	double d1_shape_2D(int node, int d, double e1, double e2);
	Segment inplanewallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k);
	Segment symplaneperiodicwallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k);

	// If a segment encounters one of the boundary planes, enforce the periodic boundary conditions
	Segment PeriodicBC(Segment &seg);

public:
	bool BC_violated;
	bool BC_bouncy;

private:
	FEAngio&	m_angio;
};
