///////////////////////////////////////////////////////////////////////
// BC.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The BC class is used by the CULTURE class to handle boundary 
// conditions within the model.
///////////////////////////////////////////////////////////////////////



#pragma once

#include "angio3d.h"
#include "Segment.h"
#include "Elem.h"
#include "Grid.h"
#include "Data.h"
#include <list>

using namespace std;

class BC
{
  public:
	BC();
	virtual ~BC();
	
	void checkBC(Segment &seg, int k, Grid &grid, Data &data, list<Segment> &frag);
	void enforceBC(vec3d i_point, int face, char bctype, Segment &seg, int elem_num, Grid &grid, int k, Data &data, list<Segment> &frag);
	void flatwallBC(vec3d i_point, int face, Segment &seg, int elem_num, int k, Grid &grid);
	Segment bouncywallBC(vec3d i_point, int face, Segment &seg, int elem_num, Grid &grid, int k, Data &data, list<Segment> &frag);
	void collfibwallBC(vec3d i_point, int face, Segment &seg, int elem_num, Grid &grid, int k, Data &data, list<Segment> &frag);
	vec3d intersceptface(int face, double &xix_0, double &xiy_0, double &xiz_0, double &xix_1, double &xiy_1, double &xiz_1, Segment &seg, int k);
	vec3d find_intersect(Elem &elem, int &face, Segment &seg, Grid &grid);
	bool search_neighbors_4_intersect(Elem &elem, int face, double &lam, double &e1, double &e2, vec3d &A, vec3d &B, Grid &grid, vec3d &inter); 
	void newton_find_intersect(double &lam, double &e1, double &e2, vec3d &A, vec3d &B, Node &X1, Node &X2, Node &X3, Node &X4);
	double shape_2D(int node, double e1, double e2);
	double d1_shape_2D(int node, int d, double e1, double e2);
	Segment inplanewallBC(vec3d i_point, int face, Segment &seg, int elem_num, Grid &grid, int k, Data &data, list<Segment> &frag);
	Segment symplaneperiodicwallBC(vec3d i_point, int face, Segment &seg, int elem_num, Grid &grid, int k, Data &data, list<Segment> &frag);

  public:
	bool BC_violated;
	bool BC_bouncy;
};
