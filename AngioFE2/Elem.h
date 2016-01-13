#pragma once
#include <FECore\vec3d.h>

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
	
	vec3d m_ecm_den_grad;
	vec3d m_collfib;		// current collagen fiber direction
	vec3d m_collfib0;		// initial collagen fiber direction

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
	bool			BC;
	unsigned int	bc_type;	// boundary condition type for this face
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
// Represents an element of the grid
class Elem
{
public:
    Elem() : elem_num(-1), volume(0.), volume0(0.), alpha(0.) {}

	// get the bounding box of this element
	BBOX GetBoundingBox();

	void GetFace(int n, int* fn);

	Node* GetNode(int i)
	{
		switch(i)
		{
		case 0: return n1; break;
		case 1: return n2; break;
		case 2: return n3; break;
		case 3: return n4; break;
		case 4: return n5; break;
		case 5: return n6; break;
		case 6: return n7; break;
		case 7: return n8; break;
		}
		return 0;
	}

	Face* GetFace(int i)
	{
		switch (i)
		{
		case 0: return &f1; break;
		case 1: return &f2; break;
		case 2: return &f3; break;
		case 3: return &f4; break;
		case 4: return &f5; break;
		case 5: return &f6; break;
		}
		return 0;
	}

public:
	int elem_num;       // Element Identifier
    
	double volume;
	double volume0;

	double alpha;
	vec3d fiber_orient;

    Node * n1;            // Bottom, lower, left node
    Node * n2;            // Bottom, lower, right node
    Node * n3;            // Bottom, upper, left node
    Node * n4;            // Bottom, upper, right node
    Node * n5;            // Top, lower, left node 
    Node * n6;            // Top, lower, right node    
    Node * n7;            // Top, upper, left node
    Node * n8;            // Top, upper, right node
    
    Face f1;            // Front face -y
    Face f2;            // Right face +x
    Face f3;            // Back face +y
    Face f4;            // Left face -x
    Face f5;            // Top face +z
    Face f6;            // Bottom face -z

	int		m_nbr[6];	// list of indices to neighbor elements (-1 if no neighbor)
};
