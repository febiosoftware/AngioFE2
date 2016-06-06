#include "stdafx.h"
#include "Grid.h"
#include "angio3d.h"
#include <FECore/FEMesh.h>
#include <FECore/vec3d.h>
#include <FECore/FESolidDomain.h>
#include <math.h>
#include "Elem.h"
#include "BC.h"

//-----------------------------------------------------------------------------
Grid::Grid(FEMesh& mesh) : m_mesh(mesh)
{
	m_coll_den	= 3.0;
	m_bc_type	= BC::STOP;
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
		node.m_id = i;													// Give the node it's ID 
		node.rt = m_mesh.Node(i).m_r0;									// Set the current position as FEBio initial position
		node.r0 = node.rt;												// Set initial position to current position

		// Add the node to the list
		m_Node.push_back(node);
	}

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
	            
			elem.m_pnode[0] = &m_Node[n1];
			elem.m_pnode[1] = &m_Node[n2];
			elem.m_pnode[2] = &m_Node[n3];
			elem.m_pnode[3] = &m_Node[n4];
			elem.m_pnode[4] = &m_Node[n5];
			elem.m_pnode[5] = &m_Node[n6];
			elem.m_pnode[6] = &m_Node[n7];
			elem.m_pnode[7] = &m_Node[n8];
	        
			m_Elem.push_back(elem);
		}
	}
	
	// Find all the element neighbors
	FindElementNeighbors();

	// build all faces
	BuildFaces();

	// initialize all element grid volumes
	InitGridVolumes();

	return true;
}

//-----------------------------------------------------------------------------
// This function is called after the node positions have been changed (i.e. after FE call)
// and updates all dependent quantities.
void Grid::Update()
{
	// loop over all nodes
	int NN = Nodes();
	for (int i = 0; i < NN; ++i)
	{
		// Update the grid node to the current position of the FE mesh
		GetNode(i).rt = m_mesh.Node(i).m_rt;
	}

	// Update the volume for each element in the grid
	// Notice that the deformation gradient is evaluated at the center only.
	int NE = Elems();
	for (int i = 0; i < NE; ++i)
	{
		Elem& elem = m_Elem[i];
		mat3d F = DeformationGradient(elem, 0.0, 0.0, 0.0);
		double Jacob = F.det();
		m_Elem[i].m_volume = Jacob*elem.m_volume0;
	}

	// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)
	update_ECM();
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
		Elem& ei = m_Elem[i];
		for (int k=0; k<6; ++k)
		{
			ei.GetFace(k, ni);

			int nbr = -1;
			for (int j=0; j<NE; ++j)
			{
				if (i != j)
				{
					Elem& ej = m_Elem[j];
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
// This function builds the face table and assigns the faces to Elem::f[i].
// This function must be called after the element neighbors are found (i.e. Grid::FindElementNeighbors())
void Grid::BuildFaces()
{
	// clear all faces
	m_Face.clear();

	// we assume that the element neighbors have been identified. 
	// we create a face for each null neighbor
	int NF = 0;
	int NE = Elems();
	for (int i=0; i<NE; ++i)
	{
		Elem& el = m_Elem[i];
		for (int j=0; j<6; ++j)
		{
			if (el.m_nbr[j] == -1) ++NF;
		}
	}

	// allocate faces
	m_Face.resize(NF);

	// assign faces
	NF = 0;
	for (int i=0; i<NE; ++i)
	{
		Elem& el = m_Elem[i];
		for (int j=0; j<6; ++j)
		{
			if (el.m_nbr[j] == -1)
			{
				Face& f = m_Face[NF];
				assert(f.m_nelem == -1);
				f.m_nelem = i;
				assert(el.m_pface[j]==0);
				el.m_pface[j] = &f;

				// assign ID
				f.id = NF++;

				// assign node numbers
				el.GetFace(j, f.m_node);

				// set the default bc type
				f.bc_type = m_bc_type;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Find the element in which the point lies.
int Grid::findelem(const vec3d& pt)
{
	const double eps = 0.00001;

	// loop over all elements
	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
		// get the next candidate
        Elem& elem = m_Elem[i];
        
 		// get the bounding box
		BBOX b = elem.GetBoundingBox();

		// inflate it a bit
		double inf = 0.01*b.Size();
		b.Inflate(inf);

		// Do a quick check to see if the point lies inside the bounding box
		if (b.IsInside(pt))
		{
			// get the natural coordinates of the point in the element
		    vec3d q(0,0,0);
			natcoord(q, pt, i);

			// see if the natural coordinates fall within the valid range
			if ((fabs(q.x) <= (1.0 + eps)) && (fabs(q.y) <= (1.0 + eps)) && (fabs(q.z) <= (1.0 + eps)))
			{
				return i;
			}
		}
	}
    return -1;
}

//-----------------------------------------------------------------------------
GridPoint Grid::FindGridPoint(int nelem, vec3d& q)
{
	GridPoint pt(nelem, q);
	pt.r = Position(pt);
	return pt;
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
    vec3d q(0,0,0);

	// loop over all elements
	int NE = Elems();
    for (int i = 0; i < NE; i++)
    {
		// get the next element
        Elem& elem = m_Elem[i];
        
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
			natcoord(q, r, i);

			// ensure that the natural coordinates lie inside the valid range [-1,1]
            if ((fabs(q.x) <= (1.0 + eps)) && (fabs(q.y) <= (1.0 + eps)) && (fabs(q.z) <= (1.0 + eps)))
			{
				p.q = q;
				p.nelem = i;
				p.r = Position(p);
				assert((p.r - r).norm() < 1e-9);
				return true;
			}
		}
	}
          
	// If we get here, we couldn't find an element
    return false;
}

//-----------------------------------------------------------------------------
// Calculate the natural coordinates from a global point (pt) and an element number.
void Grid::natcoord(vec3d& q, const vec3d& pt, int elem_num)
{
    vec3d F;
    mat3d Jmat;
    vec3d E(0,0,0);
    vec3d dE;
    
    double err = 1;
    const double tol = 1e-9;
    
    Elem& elem = m_Elem[elem_num];
    
    double dN[8][3] = {0};
    double shapeF[8] = {0};
       
    int iter = 0;
	int max_iter = 10;
	
	while ((err > tol) && (iter < max_iter)){
        double r = E.x;
        double s = E.y;
        double t = E.z;
        
        shapefunctions(shapeF, r, s, t);    
        shapefun_d1(dN, r, s ,t);

		// evaluate residual
		F = pt;
		for (int i=0; i<8; ++i)
		{
			Node& n = elem.GetNode(i);
			F.x -= shapeF[i]*n.rt.x;
			F.y -= shapeF[i]*n.rt.y;
			F.z -= shapeF[i]*n.rt.z;
		}

		// evaluate jacobian
		Jmat.zero();
		for (int i=0; i<8; ++i)
		{
			Node& n = elem.GetNode(i);
			Jmat[0][0] -= dN[i][0]*n.rt.x; Jmat[0][1] -= dN[i][1]*n.rt.x; Jmat[0][2] -= dN[i][2]*n.rt.x;
			Jmat[1][0] -= dN[i][0]*n.rt.y; Jmat[1][1] -= dN[i][1]*n.rt.y; Jmat[1][2] -= dN[i][2]*n.rt.y;
			Jmat[2][0] -= dN[i][0]*n.rt.z; Jmat[2][1] -= dN[i][1]*n.rt.z; Jmat[2][2] -= dN[i][2]*n.rt.z;
		}

		// solve linear system
        dE = Jmat.inverse()*F;
        err = dE.norm();

		// update solution
		E = E - dE;
        
		++iter;
	}
	assert(iter < max_iter);
        
    q.x = E.x;
    q.y = E.y;
    q.z = E.z;
}

//-----------------------------------------------------------------------------
// Determines the shape function values for a given position in natural coordinates
void Grid::shapefunctions(double (&shapeF)[8], double r, double s, double t)
{
    shapeF[0] = 0.125*((1 - r)*(1 - s)*(1 - t));
    shapeF[1] = 0.125*((1 + r)*(1 - s)*(1 - t));
    shapeF[2] = 0.125*((1 - r)*(1 + s)*(1 - t));
    shapeF[3] = 0.125*((1 + r)*(1 + s)*(1 - t));
    shapeF[4] = 0.125*((1 - r)*(1 - s)*(1 + t));
    shapeF[5] = 0.125*((1 + r)*(1 - s)*(1 + t));
    shapeF[6] = 0.125*((1 - r)*(1 + s)*(1 + t));
    shapeF[7] = 0.125*((1 + r)*(1 + s)*(1 + t));
}

//-----------------------------------------------------------------------------
// Calculate the shape function derivatives
void Grid::shapefun_d1(double dH[8][3], double r, double s, double t)
{
	dH[0][0] = -(1 - s)*(1 - t)*0.125; dH[0][1] = -(1 - r)*(1 - t)*0.125; dH[0][2] = -(1 - r)*(1 - s)*0.125;
	dH[1][0] =  (1 - s)*(1 - t)*0.125; dH[1][1] = -(1 + r)*(1 - t)*0.125; dH[1][2] = -(1 + r)*(1 - s)*0.125; 
	dH[2][0] = -(1 + s)*(1 - t)*0.125; dH[2][1] =  (1 - r)*(1 - t)*0.125; dH[2][2] = -(1 - r)*(1 + s)*0.125;
	dH[3][0] =  (1 + s)*(1 - t)*0.125; dH[3][1] =  (1 + r)*(1 - t)*0.125; dH[3][2] = -(1 + r)*(1 + s)*0.125;
	dH[4][0] = -(1 - s)*(1 + t)*0.125; dH[4][1] = -(1 - r)*(1 + t)*0.125; dH[4][2] =  (1 - r)*(1 - s)*0.125;
	dH[5][0] =  (1 - s)*(1 + t)*0.125; dH[5][1] = -(1 + r)*(1 + t)*0.125; dH[5][2] =  (1 + r)*(1 - s)*0.125;
	dH[6][0] = -(1 + s)*(1 + t)*0.125; dH[6][1] =  (1 - r)*(1 + t)*0.125; dH[6][2] =  (1 - r)*(1 + s)*0.125;
	dH[7][0] =  (1 + s)*(1 + t)*0.125; dH[7][1] =  (1 + r)*(1 + t)*0.125; dH[7][2] =  (1 + r)*(1 + s)*0.125;
}

//-----------------------------------------------------------------------------
// Calculate the deformation gradient tensor
mat3d Grid::DeformationGradient(Elem& elem, double r, double s, double t)
{
	// Calculate the derviative of the shape functions evaluate at each node
	double dN[8][3];
   	shapefun_d1(dN, r, s, t);

	// Calculate dX/de
	mat3d dXde; dXde.zero();
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		dXde[0][0] += n.r0.x*dN[i][0]; dXde[0][1] += n.r0.x*dN[i][1]; dXde[0][2] += n.r0.x*dN[i][2];
		dXde[1][0] += n.r0.y*dN[i][0]; dXde[1][1] += n.r0.y*dN[i][1]; dXde[1][2] += n.r0.y*dN[i][2];
		dXde[2][0] += n.r0.z*dN[i][0]; dXde[2][1] += n.r0.z*dN[i][1]; dXde[2][2] += n.r0.z*dN[i][2];
	}

	// Calculate the tensor dM, which is (dX/de)^-T * dN
	mat3d dXde_inv_trans;
	dXde_inv_trans = (dXde.inverse()).transpose();

	// Calculate F 
	mat3d F; F.zero();
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		vec3d temp(dN[i][0], dN[i][1], dN[i][2]);
		vec3d dM = dXde_inv_trans*temp;
		F[0][0] += n.rt.x*dM.x; F[0][1] += n.rt.x*dM.y; F[0][2] += n.rt.x*dM.z;
		F[1][0] += n.rt.y*dM.x; F[1][1] += n.rt.y*dM.y; F[1][2] += n.rt.y*dM.z;
		F[2][0] += n.rt.z*dM.x; F[2][1] += n.rt.z*dM.y; F[2][2] += n.rt.z*dM.z;
	}

	return F;
}

//-----------------------------------------------------------------------------
// Set the initial grid volume
// This calculates the approximate volume for each element.
void Grid::InitGridVolumes()
{
	mat3d Jacob_mat;					// Jacobian matrix
	double r = 0.0, s = 0.0, t = 0.0;	// Position of the centroid in natural coordinates
	double volume = 0.0;					// Element volume

	// Calculate the value of the derviatives of the shape functions evaulated at the centroid
	double dN[8][3];
	shapefun_d1(dN, r, s, t);

	int NE = Elems();
	for (int i = 0; i < NE; ++i){									// For each element within the mesh...
		Elem& elem = m_Elem[i];												// Obtain the element
		
		// Construct the Jacobian matrix
		Jacob_mat.zero();
		for (int j=0; j<8; ++j)
		{
			Node& n = elem.GetNode(j);
			Jacob_mat[0][0] += n.rt.x*dN[j][0]; Jacob_mat[0][1] += n.rt.x*dN[j][1]; Jacob_mat[0][2] += n.rt.x*dN[j][2];
			Jacob_mat[1][0] += n.rt.y*dN[j][0]; Jacob_mat[1][1] += n.rt.y*dN[j][1]; Jacob_mat[1][2] += n.rt.y*dN[j][2];
			Jacob_mat[2][0] += n.rt.z*dN[j][0]; Jacob_mat[2][1] += n.rt.z*dN[j][1]; Jacob_mat[2][2] += n.rt.z*dN[j][2];
		}

		if (Jacob_mat.det() != 0.)											// Calculate the volume from the Jacobian matrix using the determinant
			volume = 8*Jacob_mat.det();

		m_Elem[i].m_volume = volume;										// Set the element volume
		m_Elem[i].m_volume0 = volume;										// Set the element initial volume
	}
}

//-----------------------------------------------------------------------------
// Convert local coordinates (q) to global coordinates (pt)
vec3d Grid::nattoglobal(const vec3d& q, int elem_num)
{
    Elem& elem = m_Elem[elem_num];
    double shapeF[8];
    
    shapefunctions(shapeF, q.x, q.y, q.z);

	vec3d pt = vec3d(0,0,0);
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		pt += n.rt*shapeF[i];
	}

	return pt;
}

//-----------------------------------------------------------------------------
// Calculates the global position of a grid point.
vec3d Grid::Position(const GridPoint& pt)
{
	Elem& elem = m_Elem[pt.nelem];
    double shapeF[8];
    
    shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);

	vec3d r(0,0,0);
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		r += n.rt*shapeF[i];
	}
    
    return r;
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
bool Grid::FindIntersection(vec3d& r0, vec3d& r1, int nface, FACE_INTERSECTION& ic)
{
	// get the face
	Face& face = m_Face[nface];

	// get the face nodal coordinates
	vec3d y[4];
	y[0] = m_Node[face.m_node[0]].rt;
	y[1] = m_Node[face.m_node[1]].rt;
	y[2] = m_Node[face.m_node[2]].rt;
	y[3] = m_Node[face.m_node[3]].rt;

	// find the intersection with this face
	double lam = 0.5, r = 0.0, s = 0.0;
	if (IntersectFace(y, r0, r1, lam, r, s))
	{
		if ((lam <= 1.0001) && (lam >= 0.0))
		{
			ic.q = r0 + (r1 - r0)*(0.95*lam);
			ic.nface = nface;
			ic.nelem = face.m_nelem;
			ic.r[0] = r;
			ic.r[1] = s;
			ic.norm = FindFaceNormal(y, ic.r[0], ic.r[1]);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
// This function tries to find the face of an element that intersects the segment.
// It assumes that point1 (r0) of the segment is inside the element, whereas point2 (r1) does not.
// The intersection point is returned int q and the facet number in face.
// The function returns false if the intersection search fails.
bool Grid::FindIntersection(vec3d& r0, vec3d& r1, FACE_INTERSECTION& ic)
{
	// let's try the element's faces first
	if (ic.nelem >= 0)
	{
		Elem& el = m_Elem[ic.nelem];
		for (int i=0; i<6; ++i)
		{
			Face* pf = el.GetFace(i);
			if (pf)
			{
				if (FindIntersection(r0, r1, pf->id, ic)) return true;
			}
		}
	}

	
	// we do an expensive search over all faces.
	// TODO: set up a neighbor structure for faces 
	//       and use this to speed this up.
	int NF = Faces();
	for (int i=0; i<NF; ++i)
	{
		if (FindIntersection(r0, r1, i, ic)) return true;
	}

	assert(false);

	return false;
}

//-----------------------------------------------------------------------------
// Update the ECM field after a deformation
void Grid::update_ECM()
{
	// Natural coordinates of each node within the element
	double LUT[8][3] = {
		{-1., -1. ,-1.},
		{ 1., -1. ,-1.},
		{-1.,  1. ,-1.},
		{ 1.,  1. ,-1.},
		{-1., -1. , 1.},
		{ 1., -1. , 1.},
		{-1.,  1. , 1.},
		{ 1.,  1. , 1.}
	};

	// reset nodal data
	int NN = Nodes();
	for (int i=0; i<NN; ++i)
	{
		Node& ni = GetNode(i);
		ni.m_ntag = 0;
		ni.m_collfib = vec3d(0,0,0);
		ni.m_ecm_den = 0.0;
	}

	// For each element within the grid...
	int NE = Elems();
	for (int i = 0; i < NE; ++i)
	{
		// Obtain the element
		Elem& elem = GetElement(i);
		
		// For each node in the element...
		for (int j=0; j<8; j++)
		{
			// get the node
			Node& nj = elem.GetNode(j);

			// get the ecm density and collagen fiber
			double ecm_den  = nj.m_ecm_den0;
			vec3d coll_fib = nj.m_collfib0;
			
			// Calculate the deformation gradient tensor and jacobian at the node
			mat3d F = DeformationGradient(elem, LUT[j][0], LUT[j][1], LUT[j][2]);
			double Jacob = F.det();
			assert(Jacob > 0.0);
			
			// Update the collagen fiber orientation vector into the current configuration using F		
			coll_fib = F*coll_fib;
			coll_fib.unit();

			// Update matrix density using the Jacobian
			ecm_den = ecm_den/Jacob;

			// accumulate fiber directions and densities
			nj.m_collfib += coll_fib;
			nj.m_ecm_den += ecm_den;

			// increment counter
			nj.m_ntag++;
		}
	}

	// normalize fiber vector and average ecm density
	for (int i = 0; i < NN; ++i)
	{
		Node& ni = GetNode(i);
		assert(ni.m_ntag > 0);
		ni.m_ecm_den /= (double) ni.m_ntag;
		ni.m_collfib.unit();
	}
	
	//update_ecm_den_grad();										  // Update the ECM density gradient based on the solution from FEBio
}

//-----------------------------------------------------------------------------
// Calculates the unit direction vector of the collagen fibers at a grid point.
vec3d Grid::CollagenDirection(GridPoint& pt)
{   
	// get the element
	Elem& elem = GetElement(pt.nelem);
        
    // Obtain shape function weights
    double shapeF[8];
    shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);
    
    // Determine component of new vessel direction due to nodal collagen fiber orientation
	vec3d coll_angle(0,0,0);
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		coll_angle += n.m_collfib*shapeF[i];
	}

	// make unit vector
	coll_angle.unit();

    return coll_angle;
}

//-----------------------------------------------------------------------------
// Calculates the ecm density
double Grid::FindECMDensity(const GridPoint& pt)
{
	// get the element
	Elem& elem = GetElement(pt.nelem);
        
    // Obtain shape function weights
    double shapeF[8];
    shapefunctions(shapeF, pt.q.x, pt.q.y, pt.q.z);

	// evaluate the collagen density
	double coll_den = 0.0;
	for (int i=0; i<8; ++i)
	{
		Node& n = elem.GetNode(i);
		coll_den += n.m_ecm_den*shapeF[i];
	}

	return coll_den;
}

vec3d Grid::gradient(int elemNum, vector<double>& fn, double r, double s, double t)
{
	FESolidDomain& domain = static_cast<FESolidDomain&>(m_mesh.Domain(0));
	FESolidElement& el = domain.Element(elemNum);
	double Ji[3][3];
	domain.invjact(el, Ji, r, s, t);
	double Gr[8], Gs[8], Gt[8];
	el.shape_deriv(Gr, Gs, Gt, r, s, t);

	double Gx, Gy, Gz;

	vec3d gradf;
	int N = el.Nodes();
	for (int i=0; i<N; ++i)
	{
		// calculate global gradient of shape functions
		// note that we need the transposed of Ji, not Ji itself !
		Gx = Ji[0][0]*Gr[i]+Ji[1][0]*Gs[i]+Ji[2][0]*Gt[i];
		Gy = Ji[0][1]*Gr[i]+Ji[1][1]*Gs[i]+Ji[2][1]*Gt[i];
		Gz = Ji[0][2]*Gr[i]+Ji[1][2]*Gs[i]+Ji[2][2]*Gt[i];

		// calculate pressure gradient
		gradf.x += Gx*fn[i];
		gradf.y += Gy*fn[i];
		gradf.z += Gz*fn[i];
	}

	return gradf;
}

vec3d Grid::genericGradient(int elemNum, double Node::*material_param, double r, double s, double t)
{
	vector<double> gx = createVector(elemNum, material_param);
	return gradient(elemNum, gx, r, s, t);
}

double Grid::projectToPoint(int elemNum, vector<double>& fn, double r, double s, double t)
{
	FESolidDomain& domain = static_cast<FESolidDomain&>(m_mesh.Domain(0));
	FESolidElement& el = domain.Element(elemNum);
	double H[8];
	el.shape_fnc(H, r, s, t);
	int N = el.Nodes();
	double value = 0.0;

	for(int i=0; i<N; i++)
	{
		value += fn[i]*H[i];
	}
	return value;
}

double Grid::genericProjectToPoint(int elemNum, double Node::*material_param, double r, double s, double t)
{
	vector<double> gx = createVector(elemNum, material_param);
	return projectToPoint(elemNum, gx, r, s, t);
}

vector<double> Grid::createVector(int elemNum, double Node::*material_param)
{
	FESolidDomain& domain = static_cast<FESolidDomain&>(m_mesh.Domain(0));
	Elem& angioElem = GetElement(elemNum);
	FESolidElement& el = domain.Element(elemNum);
	int neln = el.Nodes();
	vector<double> gx(neln);

	gx[0] = angioElem.GetNode(0).*material_param;
	gx[1] = angioElem.GetNode(1).*material_param;
	gx[3] = angioElem.GetNode(2).*material_param;
	gx[2] = angioElem.GetNode(3).*material_param;
	gx[4] = angioElem.GetNode(4).*material_param;
	gx[5] = angioElem.GetNode(5).*material_param;
	gx[7] = angioElem.GetNode(6).*material_param;
	gx[6] = angioElem.GetNode(7).*material_param;

	return gx;
}