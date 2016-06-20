#pragma once
#include <FECore/vec3d.h>

//-----------------------------------------------------------------------------
class Grid;

//-----------------------------------------------------------------------------
// Represents a node on the grid
class Node
{
public:
	// default constructor
    Node();

public: 
	vec3d	r0;	// initial position of node
	vec3d	rt;	// position of node
    
	double m_ecm_den;		// current ecm density
    double m_ecm_den0;		// initial ecm density
	
	vec3d m_collfib;		// current collagen fiber direction
	vec3d m_collfib0;		// initial collagen fiber direction
	double vegf_conc;		// vegf concentration
	double m_da;
	double alpha;

	int m_id;		// node ID (zero-based)
	int m_ntag;
};

//-----------------------------------------------------------------------------
// Represents a face of the grid
class Face
{
public:
	Face();

public:
	int				id;				// face id (zero-based)
	unsigned int	bc_type;		// boundary condition type for this face
	int				m_nelem;		// element number that this face belongs to.
	int				m_node[4];		// nodal IDs
};

//-----------------------------------------------------------------------------
inline double max3(double a, double b, double c)
{
	if ((a >= b)&&(a >= c)) return a;
	if ((b >= a)&&(b >= c)) return b;
	return c;
}

//-----------------------------------------------------------------------------
struct BBOX
{
public:
	BBOX(){}
	BBOX(const vec3d& p)
	{
		xmin = xmax = p.x;
		ymin = ymax = p.y;
		zmin = zmax = p.z;
	}

	void Inflate(double e)
	{
		xmin -= e; xmax += e;
		ymin -= e; ymax += e;
		zmin -= e; zmax += e;
	}

	bool IsInside(const vec3d& p)
	{
		if ((p.x < xmin)||(p.x > xmax)) return false;
		if ((p.y < ymin)||(p.y > ymax)) return false;
		if ((p.z < zmin)||(p.z > zmax)) return false;
		return true;
	}

	void Add(const vec3d& p)
	{
		if (p.x < xmin) xmin = p.x; if (p.x > xmax) xmax = p.x;
		if (p.y < ymin) ymin = p.y; if (p.y > ymax) ymax = p.y;
		if (p.z < zmin) zmin = p.z; if (p.z > zmax) zmax = p.z;
	}

	double Size()
	{
		double wx = xmax - xmin;
		double wy = ymax - ymin;
		double wz = zmax - zmin;
		return max3(wx, wy, wz);
	}

public:
	double	xmin, xmax;
	double	ymin, ymax;
	double	zmin, zmax;
};

//-----------------------------------------------------------------------------
// Represents an element of the grid. 
// Node numbering:
// 1 = Bottom, lower, left node
// 2 = Bottom, lower, right node
// 3 = Bottom, upper, left node
// 4 = Bottom, upper, right node
// 5 = Top, lower, left node 
// 6 = Top, lower, right node    
// 7 = Top, upper, left node
// 8 = Top, upper, right node
class Elem
{
public:
    Elem();

	// get the bounding box of this element
	BBOX GetBoundingBox();

	// get the node numbers of a face
	void GetFace(int n, int* fn);

	// get a node
	Node& GetNode(int i) { return *m_pnode[i]; }

	// get a face
	Face* GetFace(int i) { return m_pface[i]; }

	// get neighbor
	int Neighbor(int i) { return m_nbr[i]; }

public:
	int elem_num;       // Element Identifier
    
	double m_volume0;	// initial (approximate) element volume
	double m_volume;	// current (approximate) element volume

	double alpha;
	vec3d fiber_orient;

private:
	Node*	m_pnode[8];		// pointer to nodes
    Face*	m_pface[6];		// pointer to boundar faces (null if face is not on the boundary)
	int		m_nbr[6];		// list of indices to neighbor elements (-1 if no neighbor)

	friend class Grid;
};
