///////////////////////////////////////////////////////////////////////
// Grid.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The GRID class stores the regular quadrilateral grid that is fit to 
// the simulation domain.  The nodes of the grid store field 
// information which influences microvessel growth.  Interpolation within
// the grid is handled by trilinear shape functions.
///////////////////////////////////////////////////////////////////////



#pragma once

#include <FECore\mat3d.h>
#include <vector>
#include "Segment.h"
#include "Elem.h"

class Data;

using namespace std;

// Define the number of nodes in the x, y, and z direction
const int xnodes = 76;          // Node distribution according to total domain size of 3822 x 2548 x 200
const int ynodes = 51;
const int znodes = 20;


class Grid  
{
public:
	///// GRID: Member Functions /////
	Grid();
	virtual ~Grid();

    virtual void create_grid();

	// GRID.IsOutsideBox - Inline function that determines if a newly created vessel Segment is outside of the domain
	bool IsOutsideBox(const Segment& seg);

	// GRID.natcoordinates - Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num);
    
	// find the element in which this point lies
    int findelem(double xpt, double ypt, double zpt);

    // GRID.shapefunctions - Determine the shape function values for a given position in natural coorindates
    void shapefunctions(double (&shapeF)[8], double xix, double xiy, double xiz);
    
    // GRID.shapefun_d1 - Determine the first order derivative of the shape functions for a particular node
    void shapefun_d1(double (&dshapeF)[3], const double xix, const double xiy, const double xiz, int node);
    
    void nattoglobal(double &xpt, double &ypt, double &zpt, double xix, double xiy, double xiz, int elem_num);

	int elem_find_neighbor(int elem_num,int neighbor_id);
	
	double find_density_scale(double coll_den);
    
public:
	// number of grid nodes
	int Nodes() { return (int) nodes.size(); }

	// number of grid elements
	int Elems() { return (int) ebin.size(); }
    
public:
	int load_cond;
	
	int xnodes;
	int ynodes;
	int znodes;
	
	double xrange[2];                                           // GRID.xrange - Array containing the minimum and maximum values in the x direction
	double yrange[2];                                           // GRID.yrange - Array containing the minimum and maximum values in the y direction
	double zrange[2];                                           // GRID.zrange - Array containing the minimum and maximum values in the z direction
	
	vector<Node> nodes;			// list of grid nodes
    vector<Elem> ebin;			// list of elements
                   
	double coll_den;
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

public:
    int iNBC;
    vector<vector<int> > ieBC;
};




   
    
    
    
///////////////////////////////////////////////////////////////////////
// Inline Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////
// IsOutsideBox
///////////////////////////////////////////////////////////////////////////////////

// GRID.IsOutsideBox - Determines if a vessel segment has grown outside of the domain
//      Input:       const Segment& seg - Newly created vessel Segment
//
//      Output:      Boolean operator 'true' if segment is outside the domain

inline bool Grid::IsOutsideBox(const Segment& seg) //determine if segment is outside box
{
	if (seg.rt[0].x < xrange[0] || seg.rt[0].x > xrange[1] ||
		seg.rt[1].x < xrange[0] || seg.rt[1].x > xrange[1] ||
		seg.rt[0].y < yrange[0] || seg.rt[0].y > yrange[1] ||
		seg.rt[1].y < yrange[0] || seg.rt[1].y > yrange[1] ||
		seg.rt[0].z < zrange[0] || seg.rt[0].z > zrange[1] ||
		seg.rt[1].z < zrange[0] || seg.rt[1].z > zrange[1])
		return true;
	else
		return false;
}
