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
// utility class for defining an intersection of a face
struct FACE_INTERSECTION 
{
	int		nelem;	// element number
	int		nface;	// face number
	double	r[2];	// natural coordinates
	vec3d	q;		// intersection point (global coordinates)
	vec3d	norm;	// noram at intersection point
};

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
	bool FindIntersection(vec3d& r0, vec3d& r1, int elem, FACE_INTERSECTION& ic);

public:
	// Inline function that determines if a newly created vessel Segment is outside of the domain
	bool IsOutsideBox(const Segment& seg);

	// Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num);

	// Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoord(vec3d& q, const vec3d& pt, int elem_num)
	{
		natcoord(q.x, q.y, q.z, pt.x, pt.y, pt.z, elem_num);
	}

	// find the element in which this point lies
	int findelem(const vec3d& pt);

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
	
	// update the grid volume
	void update_grid_volume();

	// Determine if a segment encounters one of the boundary planes, find the coordinates of the intersection point
	bool intersectPlane(Segment &Seg, int n, double intersectpt[3]);

private:
	// Find the neighbors of all elements
	void FindElementNeighbors();

	// build all the faces
	void BuildFaces();

	// initialize grid elements volume
	void InitGridVolumes();
    
public:
	// number of grid nodes
	int Nodes() { return (int) m_Node.size(); }

	// get a node
	Node& GetNode(int i) { return m_Node[i]; }

	// number of grid elements
	int Elems() { return (int) m_Elem.size(); }

	// get element
	Elem& GetElement(int i) { return m_Elem[i]; }
    
public:
	int xnodes;
	int ynodes;
	int znodes;
	
	double xrange[2];                                           // GRID.xrange - Array containing the minimum and maximum values in the x direction
	double yrange[2];                                           // GRID.yrange - Array containing the minimum and maximum values in the y direction
	double zrange[2];                                           // GRID.zrange - Array containing the minimum and maximum values in the z direction
	
	double m_coll_den;
	double den_scale;                                           // GRID.den_scale - Scaling factor based off of the collagen density within the domain
	double a, b, c;                                             // GRID.a, GRID.b, GRID.c - Parameters that describe the function that determines den_scale
	                                                            //                          given the command line input of collagen density.
	char x_bctype; 
	char y_bctype;
	char z_bctype;

	unsigned int m_bc_type;	// default bc type for faces

private:
	vector<Node> m_Node;		// list of grid nodes
    vector<Elem> m_Elem;		// list of elements
	vector<Face> m_Face;		// list of faces
                   
private:
	FEMesh&		m_mesh;
};

//-----------------------------------------------------------------------------
//      Input:       const Segment& seg - Newly created vessel Segment
//
//      Output:      Boolean operator 'true' if segment is outside the domain

inline bool Grid::IsOutsideBox(const Segment& seg) //determine if segment is outside box
{
	const vec3d& r0 = seg.tip(0).rt;
	const vec3d& r1 = seg.tip(1).rt;

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
