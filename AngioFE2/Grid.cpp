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

//-----------------------------------------------------------------------------
Grid::Grid(FEMesh& mesh) : m_mesh(mesh)
{
	m_coll_den = 3.0;

	xnodes = 0;
	ynodes = 0;
	znodes = 0;

	frontbc = 'w';
	rightbc = 'w';
	backbc = 'w';
	leftbc = 'w';
	bottombc = 'w';
	topbc = 'w';
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

	    if (ei.m_nbr[0] == -1) { ei.f1.BC = true; ei.f1.bc_type = frontbc; }
	    if (ei.m_nbr[1] == -1) { ei.f2.BC = true; ei.f2.bc_type = frontbc; }
	    if (ei.m_nbr[2] == -1) { ei.f3.BC = true; ei.f3.bc_type = frontbc; }
	    if (ei.m_nbr[3] == -1) { ei.f4.BC = true; ei.f4.bc_type = frontbc; }
	    if (ei.m_nbr[4] == -1) { ei.f5.BC = true; ei.f5.bc_type = frontbc; }
	    if (ei.m_nbr[5] == -1) { ei.f6.BC = true; ei.f6.bc_type = frontbc; }
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
int Grid::findelem(double xpt, double ypt, double zpt)
{
    int elem_num = -1;
    double bb_eps = 0.1;
	double eps = 0.00001;
    
    double xmin, xmax, ymin, ymax, zmin, zmax = {0.0};
    double xix, xiy, xiz = {0.0};

	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
        Elem& elem = ebin[i];
        
 		// get the bounding box
        xmin = elem.bb_xmin();
        xmax = elem.bb_xmax();
		ymin = elem.bb_ymin();
        ymax = elem.bb_ymax();
		zmin = elem.bb_zmin();
        zmax = elem.bb_zmax();

		// inflate it a bit
		double inf = (xmax - xmin)/100.0;
		xmin -= inf; xmax += inf;
		ymin -= inf; ymax += inf;
		zmin -= inf; zmax += inf;

        
		if ((xpt >= xmin) && (xpt <= xmax)){
            if ((ypt >= ymin) && (ypt <= ymax)){
                if ((zpt >= zmin) && (zpt <= zmax)){
                    natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                    
                    if ((fabs(xix) <= (1.0 + eps)) && (fabs(xiy) <= (1.0 + eps)) && (fabs(xiz) <= (1.0 + eps))){
                        elem_num = i;
						return elem_num;}}}}
	}

    return elem_num;
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
// TODO: Make a,b,c user settings.
// TODO: This function is only used by FEAngioMaterial. 
//       Perhaps we should move this function there.
double Grid::find_density_scale(double coll_den)
{
	if (coll_den == 3.0) return 1.0;

	// Determine the density scaling factor using the function defined by a, b, c
	double den_scale;
	double a = -0.016;
	double b = 5.1605;
	double c = 0.5112;
	den_scale = a + b*exp( -c*coll_den );

	return den_scale;
}

//-----------------------------------------------------------------------------
// This function tries to find the face of an element that intersects the segment.
// It assumes that point1 of the segment is inside the element.
// This also assumes that the facets are flat.
// TODO: generalize this for curved surfaces
bool Grid::FindIntersection(vec3d& r0, vec3d& r1, int elem, vec3d& q, int& face)
{
	vec3d d = r1 - r0;

	int nf[4];
	double lam_min = 1e99;

	face = -1;
	Elem& el = ebin[elem];
	for (int i=0; i<6; ++i)
	{
		// we only search open faces
		if (el.m_nbr[i] == -1)
		{
			el.GetFace(i, nf);
			vec3d a1 = nodes[nf[0]].rt;
			vec3d a2 = nodes[nf[1]].rt;
			vec3d a3 = nodes[nf[2]].rt;

			vec3d n = (a2-a1)^(a3-a1);

			double D = n*d;
			if (D != 0.0)
			{
				double lam = (n*(a1 - r0))/D;
				if ((lam <= 1.0001) && (lam >= 0.0) && (lam < lam_min))
				{
					q = r0 + d*(0.9*lam);
					face = i;
					lam_min = lam;
				}
			}
		}
	}

	return (face!=-1);
}
