///////////////////////////////////////////////////////////////////////
// Elem.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
//  The ELEM class contains the basic hexahedral element used to 
//  construct the grid.  The ELEM class also contains the NODE class
//  and the FACE class. 
///////////////////////////////////////////////////////////////////////



#pragma once

#include <list>
#include <FECore\vec3d.h>
#include <vector>
class Segment;
using namespace std;

//-----------------------------------------------------------------------------
// Represents a node on the grid
class Node
{
public:
	// default constructor
    Node();

	// copy constructor
	Node(const Node& n);

	// assignment operator
	void operator = (const Node& n);

	// TODO: Can we remove this?
	bool operator == (const Node& n2)
	{
		if ((rt.x == n2.rt.x) && (rt.y == n2.rt.y) && (rt.z == n2.rt.z))
			return true;
		else 
			return false;
	}

public: 
	vec3d	r0;	// initial position of node
	vec3d	rt;	// position of node
    
    double theta;
    double eta;
    
	double theta0;
	double eta0;

	double ecm_den;
    double ecm_den0;
	
	int id;
        
	bool updated;

	vec3d ecm_den_grad;
	vec3d u;
	vec3d collfib;
	vec3d collfib0;

	vector<double> ecm_den_store;
	vector<vec3d> ecm_fibril_store;
};

//-----------------------------------------------------------------------------
// Represents a face of the grid
class Face
{
public:
	Face() : BC(false), bc_type('n') {}

public:
	bool BC;
	char bc_type;
};


//-----------------------------------------------------------------------------
// Represents an element of the grid
class Elem
{
public:
    Elem() : elem_num(-1), volume(0.), volume0(0.), alpha(0.) {}
	  
	// Find the dimensions of the bounding box for the element
	// TODO: Make a class for this
    double bb_xmin();
    double bb_xmax();
    double bb_ymin();
    double bb_ymax();
	double bb_zmin();
    double bb_zmax();

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
