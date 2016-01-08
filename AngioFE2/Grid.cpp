///////////////////////////////////////////////////////////////////////
// Grid.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"

#include "Grid.h"
#include "Data.h"
#include "angio3d.h"
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
}

//-----------------------------------------------------------------------------
Grid::~Grid()
{

}

//-----------------------------------------------------------------------------
// Creates grid from FE mesh
void Grid::create_grid()
{
	int elem_num = 0;
	for (int i = 0; i < iNBC; ++i)
	{
	    elem_num = ieBC[i][0];
	        
	    if (ieBC[i][1] == 1){
	        ebin[elem_num].f1.BC = true;
	        //ebin[elem_num].f1.bc_type = y_bctype;}
			ebin[elem_num].f1.bc_type = frontbc;}
	            
	    if (ieBC[i][2] == 1){
	        ebin[elem_num].f2.BC = true;
	        //ebin[elem_num].f2.bc_type = x_bctype;}
	        ebin[elem_num].f2.bc_type = rightbc;}

	    if (ieBC[i][3] == 1){
	        ebin[elem_num].f3.BC = true;
	        //ebin[elem_num].f3.bc_type = y_bctype;}
	        ebin[elem_num].f3.bc_type = backbc;}

	    if (ieBC[i][4] == 1){
	        ebin[elem_num].f4.BC = true;
	        //ebin[elem_num].f4.bc_type = x_bctype;}
			ebin[elem_num].f4.bc_type = leftbc;}
	            
	    if (ieBC[i][5] == 1){
	        ebin[elem_num].f5.BC = true;
	        //ebin[elem_num].f5.bc_type = z_bctype;}
	        ebin[elem_num].f5.bc_type = topbc;}

	    if (ieBC[i][6] == 1){
	        ebin[elem_num].f6.BC = true;    
	        //ebin[elem_num].f6.bc_type = z_bctype;}}
			ebin[elem_num].f6.bc_type = bottombc;}
	}
}

//-----------------------------------------------------------------------------
// This is defined in FECore.
//void solve_3x3(double A[3][3], double b[3], double x[3]);

///////////////////////////////////////////////////////////////////////
// findelem
///////////////////////////////////////////////////////////////////////
/*
int Grid::findelem(double xpt, double ypt, double zpt)
{
	vector<vec3d> r(8);
	int NE = Elems();
	for (int i=0; i<NE; ++i)
	{
		// get the next element
		Elem& e = ebin[i];

		// get the element nodal coordinates
		r[0] = e.n1->rt;
		r[1] = e.n2->rt;
		r[2] = e.n4->rt;	// why are 3/4 swapped
		r[3] = e.n3->rt;
		r[4] = e.n5->rt;
		r[5] = e.n6->rt;
		r[6] = e.n8->rt;
		r[7] = e.n7->rt;	// why are 7/8 swapped

		// find the range
		double xmin = r[0].x, xmax = r[0].x;
		double ymin = r[0].y, ymax = r[0].y;
		double zmin = r[0].z, zmax = r[0].z;
		for (int j=1; j<8; ++j)
		{
			if (r[j].x < xmin) xmin = r[j].x;
			if (r[j].x > xmax) xmax = r[j].x;
			if (r[j].y < ymin) ymin = r[j].y;
			if (r[j].y > ymax) ymax = r[j].y;
			if (r[j].z < zmin) zmin = r[j].z;
			if (r[j].z > zmax) zmax = r[j].z;
		}

		// inflat a little bit
		const double inf = (xmax - xmin)/1000.0;
		xmin -= inf; xmax += inf;
		ymin -= inf; ymax += inf;
		zmin -= inf; zmax += inf;

		// see if this point is inside the box
		if ((xpt >= xmin)&&(xpt <= xmax)&&
			(ypt >= ymin)&&(ypt <= ymax)&&
			(zpt >= zmin)&&(zpt <= zmax))
		{
			// If the point y lies inside the box, we apply a Newton method to find
			// the isoparametric coordinates r
			double q[3];
			q[0] = q[1] = q[2] = 0;
			const double tol = 1e-5;
			double dr[3], norm;
			double H[8], G[8][3];
			do
			{
				H[0] = 0.125*(1 - q[0])*(1 - q[1])*(1 - q[2]);
				H[1] = 0.125*(1 + q[0])*(1 - q[1])*(1 - q[2]);
				H[2] = 0.125*(1 + q[0])*(1 + q[1])*(1 - q[2]);
				H[3] = 0.125*(1 - q[0])*(1 + q[1])*(1 - q[2]);
				H[4] = 0.125*(1 - q[0])*(1 - q[1])*(1 + q[2]);
				H[5] = 0.125*(1 + q[0])*(1 - q[1])*(1 + q[2]);
				H[6] = 0.125*(1 + q[0])*(1 + q[1])*(1 + q[2]);
				H[7] = 0.125*(1 - q[0])*(1 + q[1])*(1 + q[2]);

				G[0][0] = -0.125*(1 - q[1])*(1 - q[2]); G[0][1] = -0.125*(1 - q[0])*(1 - q[2]); G[0][2] = -0.125*(1 - q[0])*(1 - q[1]);
				G[1][0] =  0.125*(1 - q[1])*(1 - q[2]); G[1][1] = -0.125*(1 + q[0])*(1 - q[2]); G[1][2] = -0.125*(1 + q[0])*(1 - q[1]);
				G[2][0] =  0.125*(1 + q[1])*(1 - q[2]); G[2][1] =  0.125*(1 + q[0])*(1 - q[2]); G[2][2] = -0.125*(1 + q[0])*(1 + q[1]);
				G[3][0] = -0.125*(1 + q[1])*(1 - q[2]); G[3][1] =  0.125*(1 - q[0])*(1 - q[2]); G[3][2] = -0.125*(1 - q[0])*(1 + q[1]);
				G[4][0] = -0.125*(1 - q[1])*(1 + q[2]); G[4][1] = -0.125*(1 - q[0])*(1 + q[2]); G[4][2] =  0.125*(1 - q[0])*(1 - q[1]);
				G[5][0] =  0.125*(1 - q[1])*(1 + q[2]); G[5][1] = -0.125*(1 + q[0])*(1 + q[2]); G[5][2] =  0.125*(1 + q[0])*(1 - q[1]);
				G[6][0] =  0.125*(1 + q[1])*(1 + q[2]); G[6][1] =  0.125*(1 + q[0])*(1 + q[2]); G[6][2] =  0.125*(1 + q[0])*(1 + q[1]);
				G[7][0] = -0.125*(1 + q[1])*(1 + q[2]); G[7][1] =  0.125*(1 - q[0])*(1 + q[2]); G[7][2] =  0.125*(1 - q[0])*(1 + q[1]);

				double R[3] = {0}, A[3][3] = {0};
				for (int k=0; k<8; ++k)
				{
					R[0] += r[k].x*H[k];
					R[1] += r[k].y*H[k];
					R[2] += r[k].z*H[k];

					A[0][0] -= r[k].x*G[k][0]; A[0][1] -= r[k].x*G[k][1]; A[0][2] -= r[k].x*G[k][2];
					A[1][0] -= r[k].y*G[k][0]; A[1][1] -= r[k].y*G[k][1]; A[1][2] -= r[k].y*G[k][2];
					A[2][0] -= r[k].z*G[k][0]; A[2][1] -= r[k].z*G[k][1]; A[2][2] -= r[k].z*G[k][2];
				}
				R[0] = xpt - R[0];
				R[1] = ypt - R[1];
				R[2] = zpt - R[2];

				solve_3x3(A, R, dr);
				q[0] -= dr[0];
				q[1] -= dr[1];
				q[2] -= dr[2];

				norm = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
			}
			while (norm > tol);

			// see if the point r lies inside the element
			const double eps = 1.0001;
			if ((q[0] >= -eps) && (q[0] <= eps) &&
				(q[1] >= -eps) && (q[1] <= eps) && 
				(q[2] >= -eps) && (q[2] <= eps)) return i;
		}
	}
	return -1;
}
*/
int Grid::findelem(double xpt, double ypt, double zpt)
{
    int elem_num = -1;
    double bb_eps = 0.1;
	double eps = 0;
    
    double xmin, xmax, ymin, ymax, zmin, zmax = {0.0};
    double xix, xiy, xiz = {0.0};

	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
        Elem& elem = ebin[i];
        
        xmin = elem.bb_xmin()*(1. - bb_eps);
        xmax = elem.bb_xmax()*(1. + bb_eps);
		if (xmin == 0.0)
			xmin = -1.0*bb_eps*xmax;
		
		ymin = elem.bb_ymin()*(1. - bb_eps);
        ymax = elem.bb_ymax()*(1. + bb_eps);
        if (ymin == 0.0)
			ymin = -1.0*bb_eps*ymax;

		zmin = elem.bb_zmin()*(1. - bb_eps);
        zmax = elem.bb_zmax()*(1. + bb_eps);
        if (zmin == 0.0)
			zmin = -1.0*bb_eps*zmax;

        
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




double Grid::find_density_scale(double coll_den)
{
	double den_scale;
	double a = -0.016;
	double b = 5.1605;
	double c = 0.5112;
	den_scale = a + b*pow(E, -c*coll_den );                        // den_scale - Determine the density scaling factor using the function defined by a, b, c

	if (coll_den == 3.0)
		den_scale = 1.0;

	return den_scale;
}
