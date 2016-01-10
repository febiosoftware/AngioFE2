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
Grid::Grid()
{
	load_cond = 0;
	coll_den = 3.0;

	xnodes = 0;
	ynodes = 0;
	znodes = 0;

	frontbc = 'w';
	rightbc = 'w';
	backbc = 'w';
	leftbc = 'w';
	bottombc = 'w';
	topbc = 'w';

	// flatten fiber option (default to false)
	m_bzfibflat = 0;
}

//-----------------------------------------------------------------------------
Grid::~Grid()
{

}

//-----------------------------------------------------------------------------
// Creates grid from FE mesh
void Grid::CreateGrid(FEMesh &mesh, vector<vec3d>& fiber, vector<double>& density)
{
	// First, we build all the nodes
	int NN = mesh.Nodes();
	for (int i = 0; i < NN; ++i)								
	{
		Node node;														// Create a new node
		node.id = i;													// Give the node it's ID 
		node.rt = mesh.Node(i).m_r0;									// Set the current position as FEBio initial position
		node.r0 = node.rt;												// Set initial position to current position

		if(coll_den == 0.0)
		{
			node.ecm_den = density[i];									// Set the density of the ECM at the node 
			node.ecm_den0 = node.ecm_den;	
		}
		else
		{
			node.ecm_den = coll_den;									// Set the density of the ECM at the node 
			node.ecm_den0 = node.ecm_den;								// Set the initial ECM density to the current ECM density
		}

		// create a unit fiber vector denoting the local 
		// collagen orientation.
		vec3d node_fiber;
		if (load_cond == 3)
		{
			node_fiber = fiber[i];
		}
		else
		{
			node_fiber.x = 2*(float(rand())/RAND_MAX - 0.5); 
			node_fiber.y = 2*(float(rand())/RAND_MAX - 0.5); 
			node_fiber.z = 2*(float(rand())/RAND_MAX - 0.5);
		}
		
		if (m_bzfibflat == 1)
			node_fiber.z = 0.25*node_fiber.z;

		// normalize the vector
		node_fiber.unit();

		// assign the node
		node.collfib = node_fiber;

		// set the initial collagen fiber
		node.collfib0 = node.collfib;

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
	for (int d = 0; d < mesh.Domains(); d++)
	{
		FEDomain& domain = mesh.Domain(d);								// Obtain the domain from FEBio (only one domain)
		int num_elem = domain.Elements();									// Read in the total number of elements from the FEBio domain

		for (int i = 0; i < num_elem; ++i)									// Iterate through all elements...
		{
			FEElement& FEelem = domain.ElementRef(i);						// Obtain element i from the FEBio domain

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
    double tol = 1e-5;
    
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
	int max_iter = 6;
	
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



///////////////////////////////////////////////////////////////////////
// shapefunctions
///////////////////////////////////////////////////////////////////////

// GRID.shapefunctions - Determines the shape function values for a given position in natural coordinates
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



///////////////////////////////////////////////////////////////////////
// shapefun_d1
///////////////////////////////////////////////////////////////////////

// GRID.shapefun_d1 - For a given node (1 - 8), this function returns the value of the derivative of the shape function
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
	den_scale = a + b*pow(E, -c*coll_den );

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
