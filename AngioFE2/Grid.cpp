///////////////////////////////////////////////////////////////////////
// Grid.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"

#include "Grid.h"
#include "angio3d.h"
#include <FECore/FEMesh.h>
#include <FECore/vec3d.h>
#include "Elem.h"
#include "math.h"
#include "BC.h"

//-----------------------------------------------------------------------------
Grid::Grid(FEMesh& mesh) : m_mesh(mesh)
{
	m_coll_den = 3.0;

	xnodes = 0;
	ynodes = 0;
	znodes = 0;

	m_bc_type = BC::STOP;
}

//-----------------------------------------------------------------------------
Grid::~Grid()
{

}

//-----------------------------------------------------------------------------
// Creates grid from FE mesh
bool Grid::Init()
{
	// First, we build all the nodes
	int NN = m_mesh.Nodes();
	for (int i = 0; i < NN; ++i)								
	{
		Node node;														// Create a new node
		node.id = i;													// Give the node it's ID 
		node.rt = m_mesh.Node(i).m_r0;									// Set the current position as FEBio initial position
		node.r0 = node.rt;												// Set initial position to current position

		// Add the node to the list
		nodes.push_back(node);
	}

	// override the grid dimensions based on the mesh
	xrange[0] = xrange[1] = nodes[0].rt.x;
	yrange[0] = yrange[1] = nodes[0].rt.y;
	zrange[0] = zrange[1] = nodes[0].rt.z;
	for (int j=1; j<NN; ++j)
	{
		if (nodes[j].rt.x < xrange[0]) xrange[0] = nodes[j].rt.x;
		if (nodes[j].rt.x > xrange[1]) xrange[1] = nodes[j].rt.x;
		if (nodes[j].rt.y < yrange[0]) yrange[0] = nodes[j].rt.y;
		if (nodes[j].rt.y > yrange[1]) yrange[1] = nodes[j].rt.y;
		if (nodes[j].rt.z < zrange[0]) zrange[0] = nodes[j].rt.z;
		if (nodes[j].rt.z > zrange[1]) zrange[1] = nodes[j].rt.z;
	}

	// add the initial data to the stores
	for (int i = 0; i < NN; i++){
		nodes[i].ecm_den_store.push_back(nodes[i].ecm_den0);
		nodes[i].ecm_fibril_store.push_back(nodes[i].collfib0);}
	
	// Read in element connectivity from the FEBio mesh	
	for (int d = 0; d < m_mesh.Domains(); d++)
	{
		FEDomain& domain = m_mesh.Domain(d);
		int ne = domain.Elements();
		for (int i = 0; i < ne; ++i)
		{
			FEElement& FEelem = domain.ElementRef(i);

			Elem elem;
			elem.elem_num = i;
			int n1 = FEelem.m_node[0];
			int n2 = FEelem.m_node[1];
			int n3 = FEelem.m_node[3];	// Notice 2 and 3 are swapped
			int n4 = FEelem.m_node[2];
			int n5 = FEelem.m_node[4];
			int n6 = FEelem.m_node[5];
			int n7 = FEelem.m_node[7];	// Notice 6 and 7 are swapped
			int n8 = FEelem.m_node[6];
	            
			elem.n1 = &nodes[n1];
			elem.n2 = &nodes[n2];
			elem.n3 = &nodes[n3];
			elem.n4 = &nodes[n4];
			elem.n5 = &nodes[n5];
			elem.n6 = &nodes[n6];
			elem.n7 = &nodes[n7];
			elem.n8 = &nodes[n8];
	        
			ebin.push_back(elem);
		}
	}
	
	// Find all the element neighbors
	FindElementNeighbors();

	// Determine the boundary flags for the faces
	// We assume that a boundary face is determined by an element that
	// has an unassigned neighbor (i.e. -1)
	int NE = Elems();
	for (int i = 0; i < NE; ++i)
	{
	    Elem& ei = ebin[i];

	    if (ei.m_nbr[0] == -1) { ei.f1.BC = true; ei.f1.bc_type = m_bc_type; }
	    if (ei.m_nbr[1] == -1) { ei.f2.BC = true; ei.f2.bc_type = m_bc_type; }
	    if (ei.m_nbr[2] == -1) { ei.f3.BC = true; ei.f3.bc_type = m_bc_type; }
	    if (ei.m_nbr[3] == -1) { ei.f4.BC = true; ei.f4.bc_type = m_bc_type; }
	    if (ei.m_nbr[4] == -1) { ei.f5.BC = true; ei.f5.bc_type = m_bc_type; }
	    if (ei.m_nbr[5] == -1) { ei.f6.BC = true; ei.f6.bc_type = m_bc_type; }
	}

	// initialize all element grid volumes
	InitGridVolumes();

	return true;
}

//-----------------------------------------------------------------------------
// Helper function used by Grid::FindNeighbor to compare two faces.
bool face_compare(int ni[4], int nj[4])
{
	if ((ni[0]!=nj[0])&&(ni[0]!=nj[1])&&(ni[0]!=nj[2])&&(ni[0]!=nj[3])) return false;
	if ((ni[1]!=nj[0])&&(ni[1]!=nj[1])&&(ni[1]!=nj[2])&&(ni[1]!=nj[3])) return false;
	if ((ni[2]!=nj[0])&&(ni[2]!=nj[1])&&(ni[2]!=nj[2])&&(ni[2]!=nj[3])) return false;
	if ((ni[3]!=nj[0])&&(ni[3]!=nj[1])&&(ni[3]!=nj[2])&&(ni[3]!=nj[3])) return false;
	return true;
}


//-----------------------------------------------------------------------------
// This function find the neighbors of all elements
void Grid::FindElementNeighbors()
{
	int ni[4], nj[4];
	int NE = Elems();
	for (int i=0; i<NE; ++i)
	{
		Elem& ei = ebin[i];
		for (int k=0; k<6; ++k)
		{
			ei.GetFace(k, ni);

			int nbr = -1;
			for (int j=0; j<NE; ++j)
			{
				if (i != j)
				{
					Elem& ej = ebin[j];
					for (int l=0; l<6; ++l)
					{
						ej.GetFace(l, nj);
						if (face_compare(ni, nj)) { nbr = j; break; }
					}
				}

				if (nbr != -1) break;
			}

			ei.m_nbr[k] = nbr;
		}
	}
}

//-----------------------------------------------------------------------------
// Find the element in which the point lies.
int Grid::findelem(const vec3d& pt)
{
	const double eps = 0.00001;
    double xix = 0.0, xiy = 0.0, xiz = 0.0;

	// loop over all elements
	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
		// get the next candidate
        Elem& elem = ebin[i];
        
 		// get the bounding box
		BBOX b = elem.GetBoundingBox();

		// inflate it a bit
		double inf = b.Size();
		b.Inflate(inf);

		// Do a quick check to see if the point lies inside the bounding box
		if (b.IsInside(pt))
		{
			// get the natural coordinates of the point in the element
			natcoord(xix, xiy, xiz, pt.x, pt.y, pt.z, i);

			// see if the natural coordinates fall within the valid range
			if ((fabs(xix) <= (1.0 + eps)) && (fabs(xiy) <= (1.0 + eps)) && (fabs(xiz) <= (1.0 + eps)))
			{
				return i;
			}
		}
	}
    return -1;
}

//-----------------------------------------------------------------------------
// Find the GridPoint where the point pt lies. 
// This effectively finds the element number and the natural coordinates of 
// point pt inside the element.
bool Grid::FindGridPoint(const vec3d& r, GridPoint& p)
{
	// Initialize just in case the search fails
	p.nelem = -1;
	p.q = vec3d(0,0,0);

	double eps = 0.00001;
    double xix = 0.0, xiy = 0.0, xiz = 0.0;

	// loop over all elements
	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
		// get the next element
        Elem& elem = ebin[i];
        
		// Since the search for natural coordinates is expensive
		// we do a quick search to determine if this point lies
		// inside a box bounding the element.
		BBOX b = elem.GetBoundingBox();

		// inflate it a bit to overcome numerical roundoff issues
		double inf = 0.001*b.Size();
		b.Inflate(inf);
        
		// If the point lies inside the box ...
		if (b.IsInside(r))
		{
			// find the natural coordinates.
			natcoord(xix, xiy, xiz, r.x, r.y, r.z, i);

			// ensure that the natural coordinates lie inside the valid range [-1,1]
            if ((fabs(xix) <= (1.0 + eps)) && (fabs(xiy) <= (1.0 + eps)) && (fabs(xiz) <= (1.0 + eps)))
			{
				p.q = vec3d(xix, xiy, xiz);
				p.nelem = i;
				return true;
			}
		}
	}
          
	// If we get here, we couldn't find an element
    return false;
}

///////////////////////////////////////////////////////////////////////
// natcoord
///////////////////////////////////////////////////////////////////////

void Grid::natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num)
{
    vec3d F;
    mat3d Jmat;
    vec3d E;
    vec3d dE;
    vec3d newE;
    
    double err = 1;
    double tol = 1e-9;
    
    Elem elem = ebin[elem_num];
    
    double dN1[3] = {0};
    double dN2[3] = {0};
    double dN3[3] = {0};
    double dN4[3] = {0};
    double dN5[3] = {0};
    double dN6[3] = {0};
    double dN7[3] = {0};
    double dN8[3] = {0};
    double shapeF[8] = {0};
       
    int iter = 0;
	int max_iter = 10;
	
	while ((err > tol) && (iter < max_iter)){
        xix = E.x;
        xiy = E.y;
        xiz = E.z;
        
        shapefunctions(shapeF, xix, xiy, xiz);    
        shapefun_d1(dN1, xix, xiy, xiz, 1);
        shapefun_d1(dN2, xix, xiy, xiz, 2);
        shapefun_d1(dN3, xix, xiy, xiz, 3);
        shapefun_d1(dN4, xix, xiy, xiz, 4);
        shapefun_d1(dN5, xix, xiy, xiz, 5);
        shapefun_d1(dN6, xix, xiy, xiz, 6);
        shapefun_d1(dN7, xix, xiy, xiz, 7);
        shapefun_d1(dN8, xix, xiy, xiz, 8);
                
        F.x = xpt - (shapeF[0]*(*elem.n1).rt.x + shapeF[1]*(*elem.n2).rt.x + shapeF[2]*(*elem.n3).rt.x + shapeF[3]*(*elem.n4).rt.x + shapeF[4]*(*elem.n5).rt.x + shapeF[5]*(*elem.n6).rt.x + shapeF[6]*(*elem.n7).rt.x + shapeF[7]*(*elem.n8).rt.x);          
        F.y = ypt - (shapeF[0]*(*elem.n1).rt.y + shapeF[1]*(*elem.n2).rt.y + shapeF[2]*(*elem.n3).rt.y + shapeF[3]*(*elem.n4).rt.y + shapeF[4]*(*elem.n5).rt.y + shapeF[5]*(*elem.n6).rt.y + shapeF[6]*(*elem.n7).rt.y + shapeF[7]*(*elem.n8).rt.y); 
        F.z = zpt - (shapeF[0]*(*elem.n1).rt.z + shapeF[1]*(*elem.n2).rt.z + shapeF[2]*(*elem.n3).rt.z + shapeF[3]*(*elem.n4).rt.z + shapeF[4]*(*elem.n5).rt.z + shapeF[5]*(*elem.n6).rt.z + shapeF[6]*(*elem.n7).rt.z + shapeF[7]*(*elem.n8).rt.z); 
               
        Jmat[0][0] = -(dN1[0]*(*elem.n1).rt.x + dN2[0]*(*elem.n2).rt.x + dN3[0]*(*elem.n3).rt.x + dN4[0]*(*elem.n4).rt.x + dN5[0]*(*elem.n5).rt.x + dN6[0]*(*elem.n6).rt.x + dN7[0]*(*elem.n7).rt.x + dN8[0]*(*elem.n8).rt.x);
        Jmat[0][1] = -(dN1[1]*(*elem.n1).rt.x + dN2[1]*(*elem.n2).rt.x + dN3[1]*(*elem.n3).rt.x + dN4[1]*(*elem.n4).rt.x + dN5[1]*(*elem.n5).rt.x + dN6[1]*(*elem.n6).rt.x + dN7[1]*(*elem.n7).rt.x + dN8[1]*(*elem.n8).rt.x);
        Jmat[0][2] = -(dN1[2]*(*elem.n1).rt.x + dN2[2]*(*elem.n2).rt.x + dN3[2]*(*elem.n3).rt.x + dN4[2]*(*elem.n4).rt.x + dN5[2]*(*elem.n5).rt.x + dN6[2]*(*elem.n6).rt.x + dN7[2]*(*elem.n7).rt.x + dN8[2]*(*elem.n8).rt.x);
        
        Jmat[1][0] = -(dN1[0]*(*elem.n1).rt.y + dN2[0]*(*elem.n2).rt.y + dN3[0]*(*elem.n3).rt.y + dN4[0]*(*elem.n4).rt.y + dN5[0]*(*elem.n5).rt.y + dN6[0]*(*elem.n6).rt.y + dN7[0]*(*elem.n7).rt.y + dN8[0]*(*elem.n8).rt.y);
        Jmat[1][1] = -(dN1[1]*(*elem.n1).rt.y + dN2[1]*(*elem.n2).rt.y + dN3[1]*(*elem.n3).rt.y + dN4[1]*(*elem.n4).rt.y + dN5[1]*(*elem.n5).rt.y + dN6[1]*(*elem.n6).rt.y + dN7[1]*(*elem.n7).rt.y + dN8[1]*(*elem.n8).rt.y);
        Jmat[1][2] = -(dN1[2]*(*elem.n1).rt.y + dN2[2]*(*elem.n2).rt.y + dN3[2]*(*elem.n3).rt.y + dN4[2]*(*elem.n4).rt.y + dN5[2]*(*elem.n5).rt.y + dN6[2]*(*elem.n6).rt.y + dN7[2]*(*elem.n7).rt.y + dN8[2]*(*elem.n8).rt.y);
        
        Jmat[2][0] = -(dN1[0]*(*elem.n1).rt.z + dN2[0]*(*elem.n2).rt.z + dN3[0]*(*elem.n3).rt.z + dN4[0]*(*elem.n4).rt.z + dN5[0]*(*elem.n5).rt.z + dN6[0]*(*elem.n6).rt.z + dN7[0]*(*elem.n7).rt.z + dN8[0]*(*elem.n8).rt.z);
        Jmat[2][1] = -(dN1[1]*(*elem.n1).rt.z + dN2[1]*(*elem.n2).rt.z + dN3[1]*(*elem.n3).rt.z + dN4[1]*(*elem.n4).rt.z + dN5[1]*(*elem.n5).rt.z + dN6[1]*(*elem.n6).rt.z + dN7[1]*(*elem.n7).rt.z + dN8[1]*(*elem.n8).rt.z);
        Jmat[2][2] = -(dN1[2]*(*elem.n1).rt.z + dN2[2]*(*elem.n2).rt.z + dN3[2]*(*elem.n3).rt.z + dN4[2]*(*elem.n4).rt.z + dN5[2]*(*elem.n5).rt.z + dN6[2]*(*elem.n6).rt.z + dN7[2]*(*elem.n7).rt.z + dN8[2]*(*elem.n8).rt.z);
        
        dE = Jmat.inverse()*F;
        newE = E - dE;
        
        err = dE.norm();
        E = newE;
		++iter;}
        
    xix = E.x;
    xiy = E.y;
    xiz = E.z;
    
    return;
}


//-----------------------------------------------------------------------------
// Determines the shape function values for a given position in natural coordinates
//      Input:  - Position in natural coordinates (xix, xiy, xiz)
//              - Reference to array containing the 8 nodal shape function values (shapeF)
//
//      Output: - None (operates on references)

void Grid::shapefunctions(double (&shapeF)[8], double xix, double xiy, double xiz)
{
    shapeF[0] = ((1-xix)*(1-xiy)*(1-xiz))/8;                    // Shape function for node I, J, K
    shapeF[1] = ((1+xix)*(1-xiy)*(1-xiz))/8;                    // Shape function for node I+1, J, K
    shapeF[2] = ((1-xix)*(1+xiy)*(1-xiz))/8;                    // Shape function for node I, J+1, K
    shapeF[3] = ((1+xix)*(1+xiy)*(1-xiz))/8;                    // Shape function for node I+1, J+1, K
    shapeF[4] = ((1-xix)*(1-xiy)*(1+xiz))/8;                    // Shape function for node I, J, K+1
    shapeF[5] = ((1+xix)*(1-xiy)*(1+xiz))/8;                    // Shape function for node I+1, J, K+1
    shapeF[6] = ((1-xix)*(1+xiy)*(1+xiz))/8;                    // Shape function for node I, J+1, K+1
    shapeF[7] = ((1+xix)*(1+xiy)*(1+xiz))/8;                    // Shape function for node I+1, J+1, K+1
        
    return;
}

//-----------------------------------------------------------------------------
// For a given node (1 - 8), this function returns the value of the derivative of the shape function
//                    with respect to x, y, and z.
//      Input:  - Position in natural coordinates (xix, xiy, xiz)
//              - Reference to array containing the 3 first-order derivative shape function values (shapeF)
//              - Node number (1 - 8)
//
//      Output: - None (operates on references)

void Grid::shapefun_d1(double (&dshapeF)[3], const double xix, const double xiy, const double xiz, int node)
{
    if (node == 1)
    {
        dshapeF[0] = -(1 - xiy)*(1 - xiz)/8.;
        dshapeF[1] = -(1 - xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 - xix)*(1 - xiy)/8.;
    }
    
    if (node == 2)
    {
        dshapeF[0] = (1 - xiy)*(1 - xiz)/8.;
        dshapeF[1] = -(1 + xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 + xix)*(1 - xiy)/8.;
    }    
    
    if (node == 3)
    {
        dshapeF[0] = -(1 + xiy)*(1 - xiz)/8.;
        dshapeF[1] = (1 - xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 - xix)*(1 + xiy)/8.;
    }    
    
    if (node == 4)
    {
        dshapeF[0] = (1 + xiy)*(1 - xiz)/8.;
        dshapeF[1] = (1 + xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 + xix)*(1 + xiy)/8.;
    }    
    
    if (node == 5)
    {
        dshapeF[0] = -(1 - xiy)*(1 + xiz)/8.;
        dshapeF[1] = -(1 - xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 - xix)*(1 - xiy)/8.;
    }    

    if (node == 6)
    {
        dshapeF[0] = (1 - xiy)*(1 + xiz)/8.;
        dshapeF[1] = -(1 + xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 + xix)*(1 - xiy)/8.;
    }    
        
    if (node == 7)
    {
        dshapeF[0] = -(1 + xiy)*(1 + xiz)/8.;
        dshapeF[1] = (1 - xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 - xix)*(1 + xiy)/8.;
    }   

    if (node == 8)
    {
        dshapeF[0] = (1 + xiy)*(1 + xiz)/8.;
        dshapeF[1] = (1 + xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 + xix)*(1 + xiy)/8.;
    }   

    return;
}

//-----------------------------------------------------------------------------
// Calculate the deformation gradient tensor
mat3d Grid::calculate_deform_tensor(Elem& elem, const double ex, const double ey, const double ez)
{
	mat3d F;															// Deformation gradient tensor
	mat3d dXde;														// The tensor dX/de (derivative of reference position with respect to natural coordinates) 
	
	// Calculate the derviative of the shape functions evaluate at each node
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;	
   	dN1 = shapefun_d1(ex, ey, ez, 1);								
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);

	// Calculate dX/de
	dXde[0][0] = ((*elem.n1).r0.x)*dN1.x + ((*elem.n2).r0.x)*dN2.x + ((*elem.n3).r0.x)*dN3.x + ((*elem.n4).r0.x)*dN4.x + ((*elem.n5).r0.x)*dN5.x + ((*elem.n6).r0.x)*dN6.x + ((*elem.n7).r0.x)*dN7.x + ((*elem.n8).r0.x)*dN8.x;
	dXde[0][1] = ((*elem.n1).r0.x)*dN1.y + ((*elem.n2).r0.x)*dN2.y + ((*elem.n3).r0.x)*dN3.y + ((*elem.n4).r0.x)*dN4.y + ((*elem.n5).r0.x)*dN5.y + ((*elem.n6).r0.x)*dN6.y + ((*elem.n7).r0.x)*dN7.y + ((*elem.n8).r0.x)*dN8.y;
	dXde[0][2] = ((*elem.n1).r0.x)*dN1.z + ((*elem.n2).r0.x)*dN2.z + ((*elem.n3).r0.x)*dN3.z + ((*elem.n4).r0.x)*dN4.z + ((*elem.n5).r0.x)*dN5.z + ((*elem.n6).r0.x)*dN6.z + ((*elem.n7).r0.x)*dN7.z + ((*elem.n8).r0.x)*dN8.z;

	dXde[1][0] = ((*elem.n1).r0.y)*dN1.x + ((*elem.n2).r0.y)*dN2.x + ((*elem.n3).r0.y)*dN3.x + ((*elem.n4).r0.y)*dN4.x + ((*elem.n5).r0.y)*dN5.x + ((*elem.n6).r0.y)*dN6.x + ((*elem.n7).r0.y)*dN7.x + ((*elem.n8).r0.y)*dN8.x;
	dXde[1][1] = ((*elem.n1).r0.y)*dN1.y + ((*elem.n2).r0.y)*dN2.y + ((*elem.n3).r0.y)*dN3.y + ((*elem.n4).r0.y)*dN4.y + ((*elem.n5).r0.y)*dN5.y + ((*elem.n6).r0.y)*dN6.y + ((*elem.n7).r0.y)*dN7.y + ((*elem.n8).r0.y)*dN8.y;
	dXde[1][2] = ((*elem.n1).r0.y)*dN1.z + ((*elem.n2).r0.y)*dN2.z + ((*elem.n3).r0.y)*dN3.z + ((*elem.n4).r0.y)*dN4.z + ((*elem.n5).r0.y)*dN5.z + ((*elem.n6).r0.y)*dN6.z + ((*elem.n7).r0.y)*dN7.z + ((*elem.n8).r0.y)*dN8.z;

	dXde[2][0] = ((*elem.n1).r0.z)*dN1.x + ((*elem.n2).r0.z)*dN2.x + ((*elem.n3).r0.z)*dN3.x + ((*elem.n4).r0.z)*dN4.x + ((*elem.n5).r0.z)*dN5.x + ((*elem.n6).r0.z)*dN6.x + ((*elem.n7).r0.z)*dN7.x + ((*elem.n8).r0.z)*dN8.x;
	dXde[2][1] = ((*elem.n1).r0.z)*dN1.y + ((*elem.n2).r0.z)*dN2.y + ((*elem.n3).r0.z)*dN3.y + ((*elem.n4).r0.z)*dN4.y + ((*elem.n5).r0.z)*dN5.y + ((*elem.n6).r0.z)*dN6.y + ((*elem.n7).r0.z)*dN7.y + ((*elem.n8).r0.z)*dN8.y;
	dXde[2][2] = ((*elem.n1).r0.z)*dN1.z + ((*elem.n2).r0.z)*dN2.z + ((*elem.n3).r0.z)*dN3.z + ((*elem.n4).r0.z)*dN4.z + ((*elem.n5).r0.z)*dN5.z + ((*elem.n6).r0.z)*dN6.z + ((*elem.n7).r0.z)*dN7.z + ((*elem.n8).r0.z)*dN8.z;

	// Calculate the tensor dM, which is (dX/de)^-T * dN
	vec3d dM1; vec3d dM2; vec3d dM3; vec3d dM4; vec3d dM5; vec3d dM6; vec3d dM7; vec3d dM8;
	
	mat3d dXde_inv_trans;
	dXde_inv_trans = (dXde.inverse()).transpose();

	dM1 = dXde_inv_trans*dN1;
	dM2 = dXde_inv_trans*dN2;
	dM3 = dXde_inv_trans*dN3;
	dM4 = dXde_inv_trans*dN4;
	dM5 = dXde_inv_trans*dN5;
	dM6 = dXde_inv_trans*dN6;
	dM7 = dXde_inv_trans*dN7;
	dM8 = dXde_inv_trans*dN8;

	// Calculate F 
	F[0][0] = ((*elem.n1).rt.x)*dM1.x + ((*elem.n2).rt.x)*dM2.x + ((*elem.n3).rt.x)*dM3.x + ((*elem.n4).rt.x)*dM4.x + ((*elem.n5).rt.x)*dM5.x + ((*elem.n6).rt.x)*dM6.x + ((*elem.n7).rt.x)*dM7.x + ((*elem.n8).rt.x)*dM8.x;
	F[0][1] = ((*elem.n1).rt.x)*dM1.y + ((*elem.n2).rt.x)*dM2.y + ((*elem.n3).rt.x)*dM3.y + ((*elem.n4).rt.x)*dM4.y + ((*elem.n5).rt.x)*dM5.y + ((*elem.n6).rt.x)*dM6.y + ((*elem.n7).rt.x)*dM7.y + ((*elem.n8).rt.x)*dM8.y;
	F[0][2] = ((*elem.n1).rt.x)*dM1.z + ((*elem.n2).rt.x)*dM2.z + ((*elem.n3).rt.x)*dM3.z + ((*elem.n4).rt.x)*dM4.z + ((*elem.n5).rt.x)*dM5.z + ((*elem.n6).rt.x)*dM6.z + ((*elem.n7).rt.x)*dM7.z + ((*elem.n8).rt.x)*dM8.z;

	F[1][0] = ((*elem.n1).rt.y)*dM1.x + ((*elem.n2).rt.y)*dM2.x + ((*elem.n3).rt.y)*dM3.x + ((*elem.n4).rt.y)*dM4.x + ((*elem.n5).rt.y)*dM5.x + ((*elem.n6).rt.y)*dM6.x + ((*elem.n7).rt.y)*dM7.x + ((*elem.n8).rt.y)*dM8.x;
	F[1][1] = ((*elem.n1).rt.y)*dM1.y + ((*elem.n2).rt.y)*dM2.y + ((*elem.n3).rt.y)*dM3.y + ((*elem.n4).rt.y)*dM4.y + ((*elem.n5).rt.y)*dM5.y + ((*elem.n6).rt.y)*dM6.y + ((*elem.n7).rt.y)*dM7.y + ((*elem.n8).rt.y)*dM8.y;
	F[1][2] = ((*elem.n1).rt.y)*dM1.z + ((*elem.n2).rt.y)*dM2.z + ((*elem.n3).rt.y)*dM3.z + ((*elem.n4).rt.y)*dM4.z + ((*elem.n5).rt.y)*dM5.z + ((*elem.n6).rt.y)*dM6.z + ((*elem.n7).rt.y)*dM7.z + ((*elem.n8).rt.y)*dM8.z;

	F[2][0] = ((*elem.n1).rt.z)*dM1.x + ((*elem.n2).rt.z)*dM2.x + ((*elem.n3).rt.z)*dM3.x + ((*elem.n4).rt.z)*dM4.x + ((*elem.n5).rt.z)*dM5.x + ((*elem.n6).rt.z)*dM6.x + ((*elem.n7).rt.z)*dM7.x + ((*elem.n8).rt.z)*dM8.x;
	F[2][1] = ((*elem.n1).rt.z)*dM1.y + ((*elem.n2).rt.z)*dM2.y + ((*elem.n3).rt.z)*dM3.y + ((*elem.n4).rt.z)*dM4.y + ((*elem.n5).rt.z)*dM5.y + ((*elem.n6).rt.z)*dM6.y + ((*elem.n7).rt.z)*dM7.y + ((*elem.n8).rt.z)*dM8.y;
	F[2][2] = ((*elem.n1).rt.z)*dM1.z + ((*elem.n2).rt.z)*dM2.z + ((*elem.n3).rt.z)*dM3.z + ((*elem.n4).rt.z)*dM4.z + ((*elem.n5).rt.z)*dM5.z + ((*elem.n6).rt.z)*dM6.z + ((*elem.n7).rt.z)*dM7.z + ((*elem.n8).rt.z)*dM8.z;

	return F;
}

//-----------------------------------------------------------------------------
//		Evaluate the gradient of the shape functions
vec3d Grid::shapefun_d1(const double xix, const double xiy, const double xiz, int node)
{
    vec3d out;														// Output vector
	
	if (node == 1)													// dN1/de
    {
        out.x = -(1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 - xix)*(1 - xiz)/8.;
		out.z = -(1 - xix)*(1 - xiy)/8.;
    }
    
    if (node == 2)													// dN2/de
    {
        out.x = (1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 - xiy)/8.;
    }    
    
    if (node == 3)													// dN3/de
    {
        out.x = -(1 + xiy)*(1 - xiz)/8.;
        out.y = (1 - xix)*(1 - xiz)/8.;
        out.z = -(1 - xix)*(1 + xiy)/8.;
    }    
    
    if (node == 4)													// dN4/de
    {
        out.x = (1 + xiy)*(1 - xiz)/8.;
        out.y = (1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 + xiy)/8.;
    }    
    
    if (node == 5)													// dN5/de
    {
        out.x = -(1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 - xiy)/8.;
    }    

    if (node == 6)													// dN6/de
    {
        out.x = (1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 - xiy)/8.;
    }    
        
    if (node == 7)													// dN7/de
    {
        out.x = -(1 + xiy)*(1 + xiz)/8.;
        out.y = (1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 + xiy)/8.;
    }   

    if (node == 8)													// dN8/de
    {
        out.x = (1 + xiy)*(1 + xiz)/8.;
        out.y = (1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 + xiy)/8.;
    }   

    return out;
}


//-----------------------------------------------------------------------------
// Set the initial grid volume
// This calculates the approximate volume for each element.
void Grid::InitGridVolumes()
{
	mat3d Jacob_mat;														// Jacobian matrix
	double ex = 0.; double ey = 0.; double ez = 0.;						// Position of the centroid in natural coordinates
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;	// Arrays containing the derivatives of the shape functions
	double volume = 0.;													// Element volume

	// Calculate the value of the derviatives of the shape functions evaulated at the centroid
	dN1 = shapefun_d1(ex, ey, ez, 1);
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);
		
	int NE = Elems();
	for (int i = 0; i < NE; ++i){									// For each element within the mesh...
		Elem& elem = ebin[i];												// Obtain the element
		
		// Construct the Jacobian matrix
		Jacob_mat[0][0] = ((*elem.n1).rt.x)*dN1.x + ((*elem.n2).rt.x)*dN2.x + ((*elem.n3).rt.x)*dN3.x + ((*elem.n4).rt.x)*dN4.x + ((*elem.n5).rt.x)*dN5.x + ((*elem.n6).rt.x)*dN6.x + ((*elem.n7).rt.x)*dN7.x + ((*elem.n8).rt.x)*dN8.x;
		Jacob_mat[0][1] = ((*elem.n1).rt.x)*dN1.y + ((*elem.n2).rt.x)*dN2.y + ((*elem.n3).rt.x)*dN3.y + ((*elem.n4).rt.x)*dN4.y + ((*elem.n5).rt.x)*dN5.y + ((*elem.n6).rt.x)*dN6.y + ((*elem.n7).rt.x)*dN7.y + ((*elem.n8).rt.x)*dN8.y;
		Jacob_mat[0][2] = ((*elem.n1).rt.x)*dN1.z + ((*elem.n2).rt.x)*dN2.z + ((*elem.n3).rt.x)*dN3.z + ((*elem.n4).rt.x)*dN4.z + ((*elem.n5).rt.x)*dN5.z + ((*elem.n6).rt.x)*dN6.z + ((*elem.n7).rt.x)*dN7.z + ((*elem.n8).rt.x)*dN8.z;
		
		Jacob_mat[1][0] = ((*elem.n1).rt.y)*dN1.x + ((*elem.n2).rt.y)*dN2.x + ((*elem.n3).rt.y)*dN3.x + ((*elem.n4).rt.y)*dN4.x + ((*elem.n5).rt.y)*dN5.x + ((*elem.n6).rt.y)*dN6.x + ((*elem.n7).rt.y)*dN7.x + ((*elem.n8).rt.y)*dN8.x;
		Jacob_mat[1][1] = ((*elem.n1).rt.y)*dN1.y + ((*elem.n2).rt.y)*dN2.y + ((*elem.n3).rt.y)*dN3.y + ((*elem.n4).rt.y)*dN4.y + ((*elem.n5).rt.y)*dN5.y + ((*elem.n6).rt.y)*dN6.y + ((*elem.n7).rt.y)*dN7.y + ((*elem.n8).rt.y)*dN8.y;
		Jacob_mat[1][2] = ((*elem.n1).rt.y)*dN1.z + ((*elem.n2).rt.y)*dN2.z + ((*elem.n3).rt.y)*dN3.z + ((*elem.n4).rt.y)*dN4.z + ((*elem.n5).rt.y)*dN5.z + ((*elem.n6).rt.y)*dN6.z + ((*elem.n7).rt.y)*dN7.z + ((*elem.n8).rt.y)*dN8.z;
		
		Jacob_mat[2][0] = ((*elem.n1).rt.z)*dN1.x + ((*elem.n2).rt.z)*dN2.x + ((*elem.n3).rt.z)*dN3.x + ((*elem.n4).rt.z)*dN4.x + ((*elem.n5).rt.z)*dN5.x + ((*elem.n6).rt.z)*dN6.x + ((*elem.n7).rt.z)*dN7.x + ((*elem.n8).rt.z)*dN8.x;
		Jacob_mat[2][1] = ((*elem.n1).rt.z)*dN1.y + ((*elem.n2).rt.z)*dN2.y + ((*elem.n3).rt.z)*dN3.y + ((*elem.n4).rt.z)*dN4.y + ((*elem.n5).rt.z)*dN5.y + ((*elem.n6).rt.z)*dN6.y + ((*elem.n7).rt.z)*dN7.y + ((*elem.n8).rt.z)*dN8.y;
		Jacob_mat[2][2] = ((*elem.n1).rt.z)*dN1.z + ((*elem.n2).rt.z)*dN2.z + ((*elem.n3).rt.z)*dN3.z + ((*elem.n4).rt.z)*dN4.z + ((*elem.n5).rt.z)*dN5.z + ((*elem.n6).rt.z)*dN6.z + ((*elem.n7).rt.z)*dN7.z + ((*elem.n8).rt.z)*dN8.z;

		if (Jacob_mat.det() != 0.)											// Calculate the volume from the Jacobian matrix using the determinant
			volume = 8*Jacob_mat.det();

		ebin[i].volume = volume;										// Set the element volume
		ebin[i].volume0 = volume;										// Set the element initial volume
	}
}


//-----------------------------------------------------------------------------
// Update grid volume for each element after a deformation
void Grid::update_grid_volume()
{
	mat3d F;													// Deformation gradient tensor
	double Jacob = 0.;											// Jacobian (i.e., determinant of F)
	double ex = 0.; double ey = 0.; double ez = 0.;				// Position of the centroid in natural coordinates
	double new_volume = 0.;										// New element volume

	int NE = Elems();
	for (int i = 0; i < NE; ++i){								// For each element within the mesh...
		Elem& elem = ebin[i];									// Obtain the element
		F = calculate_deform_tensor(elem, ex, ey, ez);			// Calculate the deformation gradient tensor
		Jacob = F.det();										// Calculate the Jacobian by taking the determinant of F

		new_volume = Jacob*elem.volume0;						// Calculate the new element volume using the Jacobian
		ebin[i].volume = new_volume;							// Store the new element volume
	}
}


///////////////////////////////////////////////////////////////////////
// nattoglobal
///////////////////////////////////////////////////////////////////////

void Grid::nattoglobal(double &xpt, double &ypt, double &zpt, double xix, double xiy, double xiz, int elem_num)
{
    Elem elem = ebin[elem_num];
    double shapeF[8];
    
    shapefunctions(shapeF, xix, xiy, xiz);

    xpt = shapeF[0]*(*elem.n1).rt.x + shapeF[1]*(*elem.n2).rt.x + shapeF[2]*(*elem.n3).rt.x + shapeF[3]*(*elem.n4).rt.x + shapeF[4]*(*elem.n5).rt.x + shapeF[5]*(*elem.n6).rt.x + shapeF[6]*(*elem.n7).rt.x + shapeF[7]*(*elem.n8).rt.x;
    ypt = shapeF[0]*(*elem.n1).rt.y + shapeF[1]*(*elem.n2).rt.y + shapeF[2]*(*elem.n3).rt.y + shapeF[3]*(*elem.n4).rt.y + shapeF[4]*(*elem.n5).rt.y + shapeF[5]*(*elem.n6).rt.y + shapeF[6]*(*elem.n7).rt.y + shapeF[7]*(*elem.n8).rt.y;
    zpt = shapeF[0]*(*elem.n1).rt.z + shapeF[1]*(*elem.n2).rt.z + shapeF[2]*(*elem.n3).rt.z + shapeF[3]*(*elem.n4).rt.z + shapeF[4]*(*elem.n5).rt.z + shapeF[5]*(*elem.n6).rt.z + shapeF[6]*(*elem.n7).rt.z + shapeF[7]*(*elem.n8).rt.z;
    
    return;
}

//-----------------------------------------------------------------------------
// Calculates the global position of a grid point.
vec3d Grid::Position(const GridPoint& pt)
{
	Elem elem = ebin[pt.nelem];
    double shapeF[8];
    
    shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);

	vec3d r;
    r.x = shapeF[0]*(*elem.n1).rt.x + shapeF[1]*(*elem.n2).rt.x + shapeF[2]*(*elem.n3).rt.x + shapeF[3]*(*elem.n4).rt.x + shapeF[4]*(*elem.n5).rt.x + shapeF[5]*(*elem.n6).rt.x + shapeF[6]*(*elem.n7).rt.x + shapeF[7]*(*elem.n8).rt.x;
    r.y = shapeF[0]*(*elem.n1).rt.y + shapeF[1]*(*elem.n2).rt.y + shapeF[2]*(*elem.n3).rt.y + shapeF[3]*(*elem.n4).rt.y + shapeF[4]*(*elem.n5).rt.y + shapeF[5]*(*elem.n6).rt.y + shapeF[6]*(*elem.n7).rt.y + shapeF[7]*(*elem.n8).rt.y;
    r.z = shapeF[0]*(*elem.n1).rt.z + shapeF[1]*(*elem.n2).rt.z + shapeF[2]*(*elem.n3).rt.z + shapeF[3]*(*elem.n4).rt.z + shapeF[4]*(*elem.n5).rt.z + shapeF[5]*(*elem.n6).rt.z + shapeF[6]*(*elem.n7).rt.z + shapeF[7]*(*elem.n8).rt.z;
    
    return r;
}


///////////////////////////////////////////////////////////////////////
// elem_find_neighbor
///////////////////////////////////////////////////////////////////////

int Grid::elem_find_neighbor(int elem_num,int neighbor_id)
{
	int elem_neighbor = -1;
	
	if (elem_num == -1)
		return elem_neighbor;


	Elem& elem = ebin[elem_num];

	// neighbor_ids:
	// 1 - Fronty (neighbor on face with normal -y)
	// 2 - Righty (neighbor on face with normal +x)
	// 3 - Backy  (neighbor on face with noraml +y)
	// 4 - Lefty  (neighbor on face with normal -x)
	// 5 - Toppy  (neighbor on face with normal +z)
	// 6 - Bottomy (neighbor on face with normal -z)
	
	int NE = Elems();
	for (int i = 0; i < NE; i++)
	{
		if (i != elem_num)
		{
			Elem& neighbor = ebin[i];

			// Find Fronty
			if (neighbor_id == 1){
				if ((elem.n1 == neighbor.n3) && (elem.n2 == neighbor.n4) && (elem.n5 == neighbor.n7) && (elem.n6 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Righty
			if (neighbor_id == 2){
				if ((elem.n2 == neighbor.n1) && (elem.n4 == neighbor.n3) && (elem.n6 == neighbor.n5) && (elem.n8 == neighbor.n7)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Backy
			if (neighbor_id == 3){
				if ((elem.n3 == neighbor.n1) && (elem.n4 == neighbor.n2) && (elem.n7 == neighbor.n5) && (elem.n8 == neighbor.n6)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Lefty
			if (neighbor_id == 4){
				if ((elem.n1 == neighbor.n2) && (elem.n3 == neighbor.n4) && (elem.n5 == neighbor.n6) && (elem.n7 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Toppy
			if (neighbor_id == 5){
				if ((elem.n5 == neighbor.n1) && (elem.n6 == neighbor.n2) && (elem.n7 == neighbor.n3) && (elem.n8 == neighbor.n4)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Bottomy
			if (neighbor_id == 5){
				if ((elem.n1 == neighbor.n5) && (elem.n2 == neighbor.n6) && (elem.n3 == neighbor.n7) && (elem.n4 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}
		}
	}

	return elem_neighbor;
}

//-----------------------------------------------------------------------------
// calculate intersection between a line and a (quadratic) face.
// This is used by Grid::FindIntersection to determine the intersection between
// a segment and an element's face.
bool IntersectFace(vec3d p[4], const vec3d& r0, const vec3d& r1, double& lam, double& r, double& s)
{
	r = 0.0;
	s = 0.0;
	double H[4], Hr[4], Hs[4];
	vec3d R;
	mat3d K;

	vec3d dr = r1 - r0;

	const double tol = 1e-9;
	double unorm = 0.0;
	do
	{
		// evaluate shape functions
		H[0] = 0.25*(1.0 - r)*(1.0 - s);
		H[1] = 0.25*(1.0 + r)*(1.0 - s);
		H[2] = 0.25*(1.0 + r)*(1.0 + s);
		H[3] = 0.25*(1.0 - r)*(1.0 + s);

		// evaluate shape function derivatives
		Hr[0] = -0.25*(1.0 - s); Hs[0] = -0.25*(1.0 - r);
		Hr[1] =  0.25*(1.0 - s); Hs[1] = -0.25*(1.0 + r);
		Hr[2] =  0.25*(1.0 + s); Hs[2] =  0.25*(1.0 + r);
		Hr[3] = -0.25*(1.0 + s); Hs[3] =  0.25*(1.0 - r);

		vec3d yr = p[0]*Hr[0] + p[1]*Hr[1] + p[2]*Hr[2] + p[3]*Hr[3];
		vec3d ys = p[0]*Hs[0] + p[1]*Hs[1] + p[2]*Hs[2] + p[3]*Hs[3];

		// setup residual
		R = r0 + dr*lam - (p[0]*H[0] + p[1]*H[1] + p[2]*H[2] + p[3]*H[3]);

		// setup stiffness matrix
		K[0][0] = dr.x; K[0][1] = -yr.x; K[0][2] = -ys.x;
		K[1][0] = dr.y; K[1][1] = -yr.y; K[1][2] = -ys.y;
		K[2][0] = dr.z; K[2][1] = -yr.z; K[2][2] = -ys.z;

		// make sure the determinant is positive
		double D = K.det();
		if (D == 0.0) return false;

		// solve
		vec3d du = -(K.inverse()*R);
		unorm = du*du;

		// If the norm grows too big, the line will most likely not intersect with the facet
		if (unorm > 100.0) return false;

		lam += du.x;
		r += du.y;
		s += du.z; 
	}
	while (unorm > tol);

	// make sure the natural coordinates fall within range
	const double one = 1.0 + 1e-5;
	if ((r>=-one)&&(r<=one)&&(s>=-one)&&(s<=one)) return true;
	return false;
}

//-----------------------------------------------------------------------------
vec3d FindFaceNormal(vec3d y[4], double r, double s)
{
	// evaluate shape function derivatives
	double Hr[4], Hs[4];
	Hr[0] = -0.25*(1.0 - s); Hs[0] = -0.25*(1.0 - r);
	Hr[1] =  0.25*(1.0 - s); Hs[1] = -0.25*(1.0 + r);
	Hr[2] =  0.25*(1.0 + s); Hs[2] =  0.25*(1.0 + r);
	Hr[3] = -0.25*(1.0 + s); Hs[3] =  0.25*(1.0 - r);

	// calculate covariant base vectors
	vec3d dr = y[0]*Hr[0] + y[1]*Hr[1] + y[2]*Hr[2] + y[3]*Hr[3];
	vec3d ds = y[0]*Hs[0] + y[1]*Hs[1] + y[2]*Hs[2] + y[3]*Hs[3];

	// calculate unit normal
	vec3d n = dr ^ ds;
	n.unit();

	return n;
}

//-----------------------------------------------------------------------------
// This function tries to find the face of an element that intersects the segment.
// It assumes that point1 (r0) of the segment is inside the element, whereas point2 (r1) does not.
// The intersection point is returned int q and the facet number in face.
// The function returns false if the intersection search fails.
bool Grid::FindIntersection(vec3d& r0, vec3d& r1, int elem, FACE_INTERSECTION& ic)
{
	double lam_min = 1e99;

	ic.nelem = elem;
	ic.nface = -1;
	Elem& el = ebin[elem];

	vec3d d = r1 - r0;
	int nf[4];
	vec3d y[4];
	for (int i=0; i<6; ++i)
	{
		// we only search open faces
		if (el.m_nbr[i] == -1)
		{
			el.GetFace(i, nf);
			y[0] = nodes[nf[0]].rt;
			y[1] = nodes[nf[1]].rt;
			y[2] = nodes[nf[2]].rt;
			y[3] = nodes[nf[3]].rt;

			double lam = 0.5, r = 0.0, s = 0.0;
			if (IntersectFace(y, r0, r1, lam, r, s))
			{
				if ((lam <= 1.0001) && (lam >= 0.0) && (lam < lam_min))
				{
					ic.q = r0 + d*(0.95*lam);
					ic.nface = i;
					ic.r[0] = r;
					ic.r[1] = s;
					lam_min = lam;
				}
			}
		}
	}

	assert(ic.nface!=-1);

	// calculate face normal
	if (ic.nface != -1)
	{
		el.GetFace(ic.nface, nf);
		y[0] = nodes[nf[0]].rt;
		y[1] = nodes[nf[1]].rt;
		y[2] = nodes[nf[2]].rt;
		y[3] = nodes[nf[3]].rt;
		ic.norm = FindFaceNormal(y, ic.r[0], ic.r[1]);
	}

	return (ic.nface!=-1);
}

//-----------------------------------------------------------------------------
// line equation P = LP[3] + u*V[3]
// Plane equation N[3].(P[3]-P*[3]) = 0
// solving for u, u={N[3].(P*[3]-LP[3])}/{N[3].V[3]}

// box face numbering
// Front = 0
// Right = 1
// Back = 2+
// Left = 3
// Top = 4
// Bottom = 5
bool Grid::intersectPlane(Segment &seg, int n, double intersectpt[3])
{
	double N0[3], N1[3], N2[3], N3[3], N4[3], N5[3]; //face normals
	double P0[3], P1[3], P2[3], P3[3], P4[3], P5[3]; //point on faces

	//front
	N0[0] = 0;
	N0[1] = -1;
	N0[2] = 0;
	P0[0] = 0;
	P0[1] = 0;
	P0[2] = 0;

	//right
	N1[0] = 1;
	N1[1] = 0;
	N1[2] = 0;
	P1[0] = xrange[1];
	P1[1] = 0;
	P1[2] = 0;

	//back
	N2[0] = 0;
	N2[1] = 1;
	N2[2] = 0;
	P2[0] = 0;
	P2[1] = yrange[1];
	P2[2] = 0;

	//left
	N3[0] = -1;
	N3[1] = 0;
	N3[2] = 0;
	P3[0] = 0;
	P3[1] = 0;
	P3[2] = 0;

	//top
	N4[0] = 0;
	N4[1] = 0;
	N4[2] = 1;
	P4[0] = 0;
	P4[1] = 0;
	P4[2] = zrange[1];

	//bottom
	N5[0] = 0;
	N5[1] = 0;
	N5[2] = -1;
	P5[0] = 0;
	P5[1] = 0;
	P5[2] = 0;
	double V[3]; //segment displacement vector
	double V2[3]; //P*[3]-LP[3]
	double LP[3]; //origin of segment
	double u; //scalar weight to move along segment displacement vector


	if (seg.tip(1).bactive)
	{
		V[0] = seg.tip(1).rt.x - seg.tip(0).rt.x;
		V[1] = seg.tip(1).rt.y - seg.tip(0).rt.y;
		V[2] = seg.tip(1).rt.z - seg.tip(0).rt.z;
		LP[0] = seg.tip(0).rt.x;
		LP[1] = seg.tip(0).rt.y;
		LP[2] = seg.tip(0).rt.z;
	}
	else
	{
		V[0] = seg.tip(0).rt.x - seg.tip(1).rt.x;
		V[1] = seg.tip(0).rt.y - seg.tip(1).rt.y;
		V[2] = seg.tip(0).rt.z - seg.tip(1).rt.z;
		LP[0] = seg.tip(1).rt.x;
		LP[1] = seg.tip(1).rt.y;
		LP[2] = seg.tip(1).rt.z;
	}

	switch (n)
	{
		case 0: //front
			V2[0] = P0[0] - LP[0];
			V2[1] = P0[1] - LP[1];
			V2[2] = P0[2] - LP[2];
			u = vec_dot(N0,V2)/vec_dot(N0,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > xrange[0] && intersectpt[0] < xrange[1] &&
					intersectpt[2] > zrange[0] && intersectpt[2] < zrange[1])
					return true;
			}
			break;
		case 1: //right
			V2[0] = P1[0] - LP[0];
			V2[1] = P1[1] - LP[1];
			V2[2] = P1[2] - LP[2];
			u = vec_dot(N1,V2)/vec_dot(N1,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[1] > yrange[0] && intersectpt[1] < yrange[1] &&
					intersectpt[2] > zrange[0] && intersectpt[2] < zrange[1])
					return true;
			}
			break;
		case 2: //back
			V2[0] = P2[0] - LP[0];
			V2[1] = P2[1] - LP[1];
			V2[2] = P2[2] - LP[2];
			u = vec_dot(N2,V2)/vec_dot(N2,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > xrange[0] && intersectpt[0] < xrange[1] &&
					intersectpt[2] > zrange[0] && intersectpt[2] < zrange[1])
					return true;
			}
			break;
		case 3: //left
			V2[0] = P3[0] - LP[0];
			V2[1] = P3[1] - LP[1];
			V2[2] = P3[2] - LP[2];
			u = vec_dot(N3,V2)/vec_dot(N3,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[1] > yrange[0] && intersectpt[1] < yrange[1] &&
					intersectpt[2] > zrange[0] && intersectpt[2] < zrange[1])
					return true;
			}
			break;
		case 4: //top
			V2[0] = P4[0] - LP[0];
			V2[1] = P4[1] - LP[1];
			V2[2] = P4[2] - LP[2];
			u = vec_dot(N4,V2)/vec_dot(N4,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > xrange[0] && intersectpt[0] < xrange[1] &&
					intersectpt[1] > yrange[0] && intersectpt[1] < yrange[1])
					return true;
			}
			break;
		case 5: //bottom
			V2[0] = P5[0] - LP[0];
			V2[1] = P5[1] - LP[1];
			V2[2] = P5[2] - LP[2];
			u = vec_dot(N5,V2)/vec_dot(N5,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > xrange[0] && intersectpt[0] < xrange[1] &&
					intersectpt[1] > yrange[0] && intersectpt[1] < yrange[1])
					return true;
			}
			break;
	}
	return false;
}
