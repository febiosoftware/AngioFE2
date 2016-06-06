#pragma once
#include <FECore/mat3d.h>
#include <FECore/vec3d.h>
#include <FECore/FESolidDomain.h>
#include "Segment.h"
#include "Elem.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class FEMesh;

//-----------------------------------------------------------------------------
// utility class for defining an intersection of a face
struct FACE_INTERSECTION 
{
public:
	int		nelem;	// element number
	int		nface;	// face number
	double	r[2];	// natural coordinates
	vec3d	q;		// intersection point (global coordinates)
	vec3d	norm;	// noram at intersection point

public:
	FACE_INTERSECTION() { nelem = nface = -1; }
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

	// update grid (after nodes have been repositioned)
	void Update();

public:
	// number of grid nodes
	int Nodes() { return (int) m_Node.size(); }

	// get a node
	Node& GetNode(int i) { return m_Node[i]; }

	// number of grid elements
	int Elems() { return (int) m_Elem.size(); }

	// get element
	Elem& GetElement(int i) { return m_Elem[i]; }

	// number of faces
	int Faces() { return (int) m_Face.size(); }

	// get a face
	Face& GetFace(int i) { return m_Face[i]; }

public:
	// find an intersection of a segment with an element
	// The intersection point is returned in q.
	bool FindIntersection(vec3d& r0, vec3d& r1, FACE_INTERSECTION& ic);

	// Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoord(vec3d& q, const vec3d& pt, int elem_num);

	// Convert from natural to global coordinates for element
    vec3d nattoglobal(const vec3d& q, int elem_num);

	// find the element in which this point lies
	int findelem(const vec3d& pt);

	// return a GridPoint variable from a global point
	// If the point lies outside the grid, GridPoint::nelem = -1.
	bool FindGridPoint(const vec3d& r, GridPoint& p);

	// setup a grid point structure
	GridPoint FindGridPoint(int nelem, vec3d& q);

    // Determine the shape function values for a given position in natural coorindates
    void shapefunctions(double (&shapeF)[8], double r, double s, double t);
    
    // Determine the first order derivative of the shape functions for a particular node
    void shapefun_d1(double dH[8][3], double r, double s, double t);

	// Evaluate deformation gradient tensor
	mat3d DeformationGradient(Elem& elem, double r, double s, double t);

	// evaluate the spatial position of a grid point
	vec3d Position(const GridPoint& pt);

	// find the unit direction vector of the collagen
	vec3d CollagenDirection(GridPoint& pt);
	
	// Find the ecm density
	double FindECMDensity(const GridPoint& pt);

	vec3d gradient(int elemNum, vector<double>& fn, double r, double s, double t);

	vec3d genericGradient(int elemNum, double Node::*material_param, double r, double s, double t);

	double projectToPoint(int elemNum, vector<double>& fn, double r, double s, double t);

	double genericProjectToPoint(int elemNum, double Node::*material_param, double r, double s, double t);
	
private:
	// Find the neighbors of all elements
	void FindElementNeighbors();

	// build all the faces
	void BuildFaces();

	// initialize grid elements volume
	void InitGridVolumes();

	// update ECM data (after grid deformation).
	void update_ECM();

	// find the intersection with a face
	bool FindIntersection(vec3d& r0, vec3d& r1, int nface, FACE_INTERSECTION& ic);

	vector<double> createVector(int elemNum, double Node::*material_param);
    
public: // user parameters

	double m_coll_den;		// default collagen density
	double den_scale;		// Scaling factor based off of the collagen density within the domain
	double a, b, c;			// Parameters that describe the function that determines den_scale
	unsigned int m_bc_type;	// default bc type for faces

private:
	vector<Node> m_Node;		// list of grid nodes
    vector<Elem> m_Elem;		// list of elements
	vector<Face> m_Face;		// list of faces
                   
	FEMesh&		m_mesh;		// The FE mesh that generated the grid
};
