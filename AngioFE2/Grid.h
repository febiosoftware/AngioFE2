#pragma once
#include <FECore\mat3d.h>
#include "Segment.h"
#include "Elem.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class FEMesh;

//-----------------------------------------------------------------------------
// Define the number of nodes in the x, y, and z direction
const int xnodes = 76;          // Node distribution according to total domain size of 3822 x 2548 x 200
const int ynodes = 51;
const int znodes = 20;

//-----------------------------------------------------------------------------
// The GRID class stores the regular quadrilateral grid that is fit to 
// the simulation domain.  The nodes of the grid store field 
// information which influences microvessel growth.  Interpolation within
// the grid is handled by trilinear shape functions.
class Grid  
{
public:
	Grid(FEMesh& mesh);
	virtual ~Grid();

	// Create the grid from the FEMesh
    bool Init();

	// find an intersection of a segment with an element
	// The intersection point is returned in q.
	bool FindIntersection(vec3d& r0, vec3d& r1, int elem, vec3d& q, int& face);

public:
	// Inline function that determines if a newly created vessel Segment is outside of the domain
	bool IsOutsideBox(const Segment& seg);

	// Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num);
    
	// find the element in which this point lies
    int findelem(double xpt, double ypt, double zpt);

	// find the element in which this point lies
	int findelem(vec3d& pt) { return findelem(pt.x, pt.y, pt.z); }

	// return a GridPoint variable from a global point
	// If the point lies outside the grid, GridPoint::nelem = -1.
	bool FindGridPoint(const vec3d& r, GridPoint& p);

    // Determine the shape function values for a given position in natural coorindates
    void shapefunctions(double (&shapeF)[8], double xix, double xiy, double xiz);
    
    // Determine the first order derivative of the shape functions for a particular node
    void shapefun_d1(double (&dshapeF)[3], const double xix, const double xiy, const double xiz, int node);

	// Another function for evaluation shape function derivatives
	// TODO: I think this function basically does the same as the one above. Perhaps one can be eliminated.
	vec3d shapefun_d1(const double xix, const double xiy, const double xiz, int node);

	// Evaluate deformation gradient tensor
	mat3d calculate_deform_tensor(Elem& elem, double e1, double e2, double e3);
    
	// Convert from natural to global coordinates for element
    void nattoglobal(double &xpt, double &ypt, double &zpt, double xix, double xiy, double xiz, int elem_num);

	vec3d Position(const GridPoint& pt);

	// find the neighbor of an element
	int elem_find_neighbor(int elem_num,int neighbor_id);
	
	// find the density scale
	double find_density_scale(double coll_den);

	// update the grid volume
	void update_grid_volume();

private:
	// Find the neighbors of all elements
	void FindElementNeighbors();

	// initialize grid elements volume
	void InitGridVolumes();
    
public:
	// number of grid nodes
	int Nodes() { return (int) nodes.size(); }

	// number of grid elements
	int Elems() { return (int) ebin.size(); }
    
public:
	int xnodes;
	int ynodes;
	int znodes;
	
	double xrange[2];                                           // GRID.xrange - Array containing the minimum and maximum values in the x direction
	double yrange[2];                                           // GRID.yrange - Array containing the minimum and maximum values in the y direction
	double zrange[2];                                           // GRID.zrange - Array containing the minimum and maximum values in the z direction
	
	vector<Node> nodes;			// list of grid nodes
    vector<Elem> ebin;			// list of elements
                   
	double m_coll_den;
	double den_scale;                                           // GRID.den_scale - Scaling factor based off of the collagen density within the domain
	double a, b, c;                                             // GRID.a, GRID.b, GRID.c - Parameters that describe the function that determines den_scale
	                                                            //                          given the command line input of collagen density.
	char x_bctype; 
	char y_bctype;
	char z_bctype;

	char frontbc;
	char rightbc;
	char backbc;
	char leftbc;
	char bottombc;
	char topbc;

private:
	FEMesh&		m_mesh;
};

//-----------------------------------------------------------------------------
//      Input:       const Segment& seg - Newly created vessel Segment
//
//      Output:      Boolean operator 'true' if segment is outside the domain

inline bool Grid::IsOutsideBox(const Segment& seg) //determine if segment is outside box
{
	const vec3d& r0 = seg.m_tip[0].rt;
	const vec3d& r1 = seg.m_tip[1].rt;

	if (r0.x < xrange[0] || r0.x > xrange[1] ||
		r1.x < xrange[0] || r1.x > xrange[1] ||
		r0.y < yrange[0] || r0.y > yrange[1] ||
		r1.y < yrange[0] || r1.y > yrange[1] ||
		r0.z < zrange[0] || r0.z > zrange[1] ||
		r1.z < zrange[0] || r1.z > zrange[1])
		return true;
	else
		return false;
}
