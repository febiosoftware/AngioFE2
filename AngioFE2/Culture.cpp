#include "stdafx.h"
#include <iostream>
#include "Culture.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include "angio3d.h"

//-----------------------------------------------------------------------------
Culture::Culture(FEAngio& angio) : bc(angio), m_angio(angio)
{
	W[0] = 0;
	W[1] = 0;
	W[2] = 0;
	W[3] = 0;

	NFRAGS = 3;
    m_num_vessel = NFRAGS - 1;

	m_num_zdead = 0;

	// Initialize anastomoses counter
	m_num_anastom = 0;
}

//-----------------------------------------------------------------------------
Culture::~Culture()
{

}

//-----------------------------------------------------------------------------
// Create initial fragments
void Culture::SeedFragments()
{
	for (int i=0; i < NFRAGS; ++i)
	{
		// Create an initial segment
		Segment seg = createInitFrag();

		// Give the segment the appropriate label
		seg.label = i;

		// Set the segment vessel as the segment label
		seg.vessel = seg.label;

		// Store the segment's time of birth
		seg.TofBirth = m_angio.m_t;

		// increment the number of segments
		m_angio.m_nsegs += 1;
		seg.seg_num = m_angio.m_nsegs;

		// add it to the list
		m_frag.push_back(seg);
	}
}

//-----------------------------------------------------------------------------
//      Input:  - DATA object
//              - GRID object
//              - Segment container (frag)
//
//      Output: - Initial fragment (seg)
//
Segment Culture::createInitFrag()
{
	Grid& grid = m_angio.grid;

    double xix, xiy, xiz = {0.};
    double xpt, ypt, zpt = {0.};

	// Create a segment
    Segment seg;

	// Set seg length to value of growth function at t = 0
	seg.length = m_angio.m_d;
    
	int elem_num = -1;
	do
	{
		// --- position randomly
		// pick an element
		float f = float(rand())/RAND_MAX;
		elem_num = int(f*grid.Elems());
		assert((elem_num>=0)&&(elem_num<grid.Elems()));

		// generate random natural coordinates
		xix = 2*((float(rand())/RAND_MAX) - 0.5);
		xiy = 2*((float(rand())/RAND_MAX) - 0.5);
		xiz = 2*((float(rand())/RAND_MAX) - 0.5);

		// convert to global coordinates
		grid.nattoglobal(xpt, ypt, zpt, xix, xiy, xiz, elem_num);
	
		// set the position of the first tip
		seg.m_tip[0].rt = vec3d(xpt, ypt, zpt);
		seg.m_tip[0].elem = elem_num;
    
		// Determine vessel orientation based off of collagen fiber orientation
		seg.uvect = findCollAngle(xpt, ypt, zpt);

		// End of new segment is origin plus length component in each direction	
		seg.m_tip[1].rt = seg.m_tip[0].rt + seg.uvect*seg.length;				  // Determine the x-coordinate of the end point using the length vector and orientation angles                 

		// make the tips active
		seg.m_tip[0].active = -1;                                            // Set the tip at the start point of the segment as -1 
		seg.m_tip[1].active = 1;                                             // Set the tip at the end point of the segment as +1
	
		// Set sprout as an initial fragment
		seg.m_sprout = Segment::SPROUT_INIT;
		
		// find the element where the second tip is
		elem_num = grid.findelem(seg.m_tip[1].rt);
	}
	while (elem_num == -1);

	seg.m_tip[1].elem = elem_num;

	return seg;
}



///////////////////////////////////////////////////////////////////////
// createNewSeg
///////////////////////////////////////////////////////////////////////

// CULTURE.createNewSeg - Create a new segment at the tip of an existing segment
//      Input:  - Iterator for Segment container 'frag', points to the parent segment (it)
//              - GRID object
//              - DATA object
//              - Index indicating which tip of the parent segment is forming the new segment (k)
//              - Segment container (frag)
//
//      Output: - Newly created segment (seg)

Segment Culture::createNewSeg(list<Segment>::iterator it, int k)
{
	Grid& grid = m_angio.grid;

	Segment seg;                                                // Declare SEGMENT seg
	seg.Recent_branch = it->Recent_branch/2;                    // Halve the Recent_branch indicator
    
	++m_angio.m_nsegs;                                               // Iterate the total segment counter +1 
	seg.seg_num = m_angio.m_nsegs;				
	
	//seg.line_num = it->line_num;

    int elem_num = 0;
    
	if (it->m_tip[k].active == 1)                                          // If the parent vessel active tip is set as +1...
	{
		//seg.length = findLength(it->x[1],it->y[1],it->z[1],grid,data);      // Determine length of new segment                   
		seg.length = m_angio.m_vess_length;

	
		double den_scale = 1.0;
		den_scale = findDenScale(it->m_tip[1].rt.x, it->m_tip[1].rt.y, it->m_tip[1].rt.z);
			
		//seg.phi1 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,1);    // Determine the angle between the new segment and the x-axis
		//seg.phi2 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,2);    // Determine the angle between the new segment and the x-y plane
		seg.uvect = findAngle(it,it->m_tip[1].rt.x,it->m_tip[1].rt.y,it->m_tip[1].rt.z); 
		
		seg.length = den_scale*m_angio.m_vess_length;

		//if (data.branch)                                                   // If new segment is a branch...
		//{
		//	if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//			seg.phi1 = ang_dom(it->phi1+pi/2);
		//        else
		//	        seg.phi1 = ang_dom(it->phi1-pi/2);
		//	
		//	//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//	//		seg.phi2 = ang_dom(it->phi2+pi/2);
		// //       else
		//	//        seg.phi2 = ang_dom(it->phi2-pi/2);
		//	
		//	//seg.phi2 = 0.1*seg.phi2;                                           // force the branch to be mostly planar (x-y plane)
		//}
		//
		
		if (m_angio.m_branch)                                                   // If new segment is a branch...
		{
			//double x = seg.uvect.x;
			//double y = seg.uvect.y;
			//double H = sqrt(x*x + y*y);
			//double theta = acos(x/H);

			//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
			//		theta = theta + pi/2;
		 //       else
			//        theta = theta - pi/2;
			//
			//seg.uvect.x = H*cos(theta);
			//seg.uvect.y = H*sin(theta);
			////seg.uvect.z = 0;

			//seg.uvect = seg.uvect/seg.uvect.norm();
		
			//vec3d coll_fib = findCollAngle(it->x[1], it->y[1], it->z[1], grid, data);

			//if (acos(coll_fib*seg.uvect) > pi/2)
			//	coll_fib = -coll_fib;

			//seg.uvect = (seg.uvect + coll_fib)/2;
			
			vec3d coll_fib = findCollAngle(it->m_tip[1].rt.x, it->m_tip[1].rt.y, it->m_tip[1].rt.z);
			vec3d newseg;
			newseg = coll_fib - seg.uvect*(seg.uvect*coll_fib)*0.5;
			newseg = newseg/newseg.norm();
			seg.uvect = newseg;
		}

		seg.label = it->label;                                  // Transfer the label of the parent segment to the new segment
		seg.vessel = it->vessel;                                // Transfer the vessel number of the parent segment to the new segment	
		
		seg.m_tip[0].rt = it->m_tip[1].rt;                                    // Set the origin of new segment as the active tip of the previous segment
		seg.m_tip[0].elem = it->m_tip[1].elem;
		seg.seg_conn[0][0] = it->seg_num;

		//seg.x[1] = seg.x[0]+seg.length*cos(seg.phi2)*cos(seg.phi1);     // Determine the x-coordinate of the end point using the length vector and orientation angles
		//seg.y[1] = seg.y[0]+seg.length*cos(seg.phi2)*sin(seg.phi1);     // Determine the y-coordinate of the end point using the length vector and orientation angles
		//seg.z[1] = seg.z[0]+seg.length*sin(seg.phi2);                   // Determine the z-coordinate of the end point using the length vector and orientation angles
    
		seg.m_tip[1].rt = seg.m_tip[0].rt + seg.uvect*seg.length;					  // Determine the x-coordinate of the end point using the length vector and orientation angles
		seg.findlength();

        seg.m_tip[1].active = 1;                                         // Turn on end tip of new segment
		seg.m_tip[0].active = 0;                                         // Turn off origin tip of new segment
		
		seg.TofBirth = m_angio.m_t;                                  // Stamp segment with time of birth
		it->m_tip[k].active = 0;                                         // Turn off previous segment tip
		seg.m_tip[k].bdyf_id = it->m_tip[k].bdyf_id;
        
		seg.m_sprout = Segment::SPROUT_POS;                                         // Set sprout for the new segment as +1, indicating this segment originated from a +1 tip
					    
		if (it->seg_conn[1][0] == 0)
			it->seg_conn[1][0] = seg.seg_num;
		else
			it->seg_conn[1][1] = seg.seg_num;
				
		elem_num = grid.findelem(seg.m_tip[1].rt);
		
		if (elem_num == -1)
		{
			bc.checkBC(seg, 1);
		}
		else
			seg.m_tip[1].elem = elem_num;

	}
	
	
	else if (it->m_tip[k].active == -1)                                    // If the parent vessel active tip is set as +1...
	{
	    //seg.length = -1*findLength(it->x[0],it->y[0],it->z[0],grid,data);  // Determine length of new segment 
		seg.length = -m_angio.m_vess_length;
		
		double den_scale = 1.0;
		den_scale = findDenScale(it->m_tip[0].rt.x, it->m_tip[0].rt.y, it->m_tip[0].rt.z);
			
		//seg.phi1 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,1);    // Determine the angle between the new segment and the x-axis
		//seg.phi2 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,2);    // Determine the angle between the new segment and the x-y plane
		
		seg.uvect = findAngle(it,it->m_tip[0].rt.x,it->m_tip[0].rt.y,it->m_tip[0].rt.z);

		seg.length = -den_scale*m_angio.m_vess_length;
				
		if (m_angio.m_branch)                                                   // If new segment is a branch...
		{
			vec3d coll_fib = findCollAngle(it->m_tip[1].rt.x, it->m_tip[1].rt.y, it->m_tip[1].rt.z);
			vec3d newseg;
			newseg = coll_fib - seg.uvect*(seg.uvect*coll_fib)*0.5;
			newseg = newseg/newseg.norm();

			seg.uvect = newseg;
		}

		seg.label = it->label;                                  // Transfer the label of the parent segment to the new segment
		seg.vessel = it->vessel;                                // Transfer the vessel number of the parent segment to the new segment
		
		seg.m_tip[1].rt = it->m_tip[0].rt;                                    // Set the origin of new segment as the active tip of the previous segment
		seg.m_tip[1].elem = it->m_tip[0].elem;
		seg.seg_conn[1][0] = it->seg_num;
		
		//seg.x[0] = seg.x[1]+seg.length*cos(seg.phi2)*cos(seg.phi1);     // Determine the x-coordinate of the end point using the length vector and orientation angles
		//seg.y[0] = seg.y[1]+seg.length*cos(seg.phi2)*sin(seg.phi1);     // Determine the y-coordinate of the end point using the length vector and orientation angles
		//seg.z[0] = seg.z[1]+seg.length*sin(seg.phi2);                   // Determine the z-coordinate of the end point using the length vector and orientation angles
		
		seg.m_tip[0].rt = seg.m_tip[1].rt + seg.uvect*seg.length;					  // Determine the x-coordinate of the end point using the length vector and orientation angles
		seg.findlength();

		seg.m_tip[0].active = -1;                                        // Turn on end tip of new segment
		seg.m_tip[1].active = 0;                                         // Turn off origin tip of new segment
		
		seg.TofBirth = m_angio.m_t;                                  // Stamp segment with time of birth
		it->m_tip[k].active = 0;                                         // Turn off previous segment tip
		seg.m_tip[k].bdyf_id = it->m_tip[k].bdyf_id;

		seg.m_sprout = Segment::SPROUT_NEG;                                        // Set sprout for the new segment as -1, indicating this segment originated from a -1 tip

		if (it->seg_conn[0][0] == 0)
			it->seg_conn[0][0] = seg.seg_num;
		else
			it->seg_conn[0][1] = seg.seg_num;
		
		
		elem_num = grid.findelem(seg.m_tip[0].rt);
		
		if (elem_num == -1)
			bc.checkBC(seg, 0);
		else
			seg.m_tip[0].elem = elem_num;
	}
	
	return seg;                                                 // Return the new segment 
}

///////////////////////////////////////////////////////////////////////
// findAngle
///////////////////////////////////////////////////////////////////////

// CULTURE.findAngle - Determine the orientation angle of a newly created segment
//      Input:  - Iterator for segment container which points to the parent segment (it)
//              - Coordinates of the active tip that is sprouting the new segment (xpt, ypt, zpt)
//              - GRID object
//              - DATA object
//              - Integer describing the angle type: 1 for phi1, 2 for phi2 (ang_type)
//
//      Output: - Segment orientation angle phi1 or phi1 (in radians)

vec3d Culture::findAngle(list<Segment>::iterator it, double xpt, double ypt, double zpt)
{
	vec3d angle;

	double den_scale = 1.0;
	
	den_scale = findDenScale(xpt, ypt, zpt);

	//double W[4] = {10*(1/grid.den_scale), 0, 0, 100};
	//double W[4] = {10, 0, 0, 100};                       // W[0] = Weight for collagen orientation
	//double W[4] = {10, 0, 0, 50};                                                                        // W[1] = Weight for vessel density
	                                                                        // W[2] = Weight for random component
	                                                                        // W[3] = Weight for previous vessel direction
	                                        
    
    vec3d coll_angle;      // Component of new vessel orientation resulting from collagen fiber orientation
    vec3d den_angle;       // Component of new vessel orientation resulting from vessel density
    vec3d ran_angle;       // Component of new vessel orientation resulting from random walk
    vec3d per_angle;       // Component of new vessel orientation resulting from previous vessel direction        
        
    
    // Find the component of the new vessel direction determined by collagen fiber orientation    
    coll_angle = findCollAngle(xpt, ypt, zpt);

	per_angle = it->uvect;

	angle = (coll_angle*W[0] + per_angle*W[3])/(W[0]+W[3]);

	angle = angle/angle.norm();

	return angle;
	
}



///////////////////////////////////////////////////////////////////////
// findCollAngle
///////////////////////////////////////////////////////////////////////


vec3d Culture::findCollAngle(double xpt, double ypt, double zpt)
{   
	Grid& grid = m_angio.grid;
    
    vec3d coll_angle;

    double xix, xiy, xiz;
    double shapeF[8];
    
    int elem_num = grid.findelem(xpt, ypt, zpt);
    
    if (elem_num < 0)
		return coll_angle;
		
	Elem elem;
    elem = grid.ebin[elem_num];
        
    // Convert to natural coordinates -1 <= Xi <= +1
    grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);        
    
    // Obtain shape function weights
    grid.shapefunctions(shapeF, xix, xiy, xiz);
    
    // Determine component of new vessel direction due to nodal collagen fiber orientation
	coll_angle = ((*elem.n1).collfib)*shapeF[0] + ((*elem.n2).collfib)*shapeF[1] + ((*elem.n3).collfib)*shapeF[2] + ((*elem.n4).collfib)*shapeF[3] + ((*elem.n5).collfib)*shapeF[4] + ((*elem.n6).collfib)*shapeF[5] + ((*elem.n7).collfib)*shapeF[6] + ((*elem.n8).collfib)*shapeF[7];
	
	coll_angle = coll_angle/coll_angle.norm();

    return coll_angle;
}



///////////////////////////////////////////////////////////////////////
// findDenScale
///////////////////////////////////////////////////////////////////////

double Culture::findDenScale(double xpt, double ypt, double zpt)
{
	Grid& grid = m_angio.grid;

	double coll_den = 0.0;
    double den_scale = 1.0;

    double xix, xiy, xiz;
    double shapeF[8];

	int elem_num = grid.findelem(xpt, ypt, zpt);
    
    if (elem_num < 0)
		return den_scale;
		
	Elem elem;
    elem = grid.ebin[elem_num];
        
    // Convert to natural coordinates -1 <= Xi <= +1
    grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);        
    
    // Obtain shape function weights
    grid.shapefunctions(shapeF, xix, xiy, xiz);
    
	coll_den = shapeF[0]*(*elem.n1).ecm_den + shapeF[1]*(*elem.n2).ecm_den + shapeF[2]*(*elem.n3).ecm_den + shapeF[3]*(*elem.n4).ecm_den + shapeF[4]*(*elem.n5).ecm_den + shapeF[5]*(*elem.n6).ecm_den + shapeF[6]*(*elem.n7).ecm_den + shapeF[7]*(*elem.n8).ecm_den;
    
	den_scale = grid.find_density_scale(coll_den);

    if (den_scale < 0.)
		den_scale = 0.;
	
	return den_scale;
}



///////////////////////////////////////////////////////////////////////
// connectSegment
///////////////////////////////////////////////////////////////////////

// creates a new segment to connect close segments

Segment Culture::connectSegment(list<Segment>::iterator it,list<Segment>::iterator it2, int k, int kk)
 {
 	Segment seg;
 	
	seg.length = (it->m_tip[k].rt - it2->m_tip[kk].rt).norm();
	
 	seg.m_tip[0].rt = it->m_tip[k].rt;
	seg.m_tip[0].elem = it->m_tip[k].elem;
	seg.m_tip[1].rt = it2->m_tip[kk].rt;
	seg.m_tip[1].elem = it2->m_tip[kk].elem;
 	
	seg.TofBirth = m_angio.m_t;
 	seg.label = it->label;
 	seg.vessel = it->vessel;
 	
 	seg.m_tip[0].active = 0;
 	seg.m_tip[1].active = 0;
 	seg.anast = 1;
 	it->anast = 1;
 	it2->anast = 1;
	seg.findlength();

 	++m_num_anastom;
 	
	++m_angio.m_nsegs;
	seg.seg_num = m_angio.m_nsegs;

	seg.seg_conn[0][0] = it->seg_num;
	
	if (k == 0)
		it->seg_conn[0][0] = seg.seg_num;
	else
		it->seg_conn[1][0] = seg.seg_num;

	seg.seg_conn[1][0] = it2->seg_num;

	if (kk == 0)
		it2->seg_conn[0][1] = seg.seg_num;
	else
		it2->seg_conn[1][1] = seg.seg_num;
	
	return seg;
 }
 
///////////////////////////////////////////////////////////////////////
// CheckForIntersection
///////////////////////////////////////////////////////////////////////
// Description: checks for intersection between a passed segment and all other existing segments that are not members
// of vessel containing the segment

void Culture::CheckForIntersection(Segment &seg, list<Segment>::iterator it)
{
	double p1[3], p2[3], pp1[3], pp2[3]; //tip points for segment to check intersection
	list<Segment>::iterator it2; //loop over segments in list
	double intersectpt[3]; //the intersection pt of vectors
	
	intersectpt[0] = 0;
	intersectpt[1] = 0;
	
	double lambda;

	vec3d& r0 = seg.m_tip[0].rt;
	vec3d& r1 = seg.m_tip[1].rt;
	
	if (seg.m_tip[1].active == 1)
	{
		p1[0] = r0.x;
		p2[0] = r1.x;
		p1[1] = r0.y;
		p2[1] = r1.y;
		p1[2] = r0.z;
		p2[2] = r1.z;
	}
	else if (seg.m_tip[0].active == -1)
	{
		p1[0] = r1.x;
		p2[0] = r0.x;
		p1[1] = r1.y;
		p2[1] = r0.y;
		p1[2] = r1.z;
		p2[2] = r0.z;
	}

	for (it2 = m_angio.cult.m_frag.begin(); it2 != m_angio.cult.m_frag.end(); ++it2)
	{
		if (it->label!=it2->label)
		{
			pp1[0] = r0.x;
			pp1[1] = r0.y;
			pp1[2] = r0.z;
			pp2[0] = r1.x;
			pp2[1] = r1.y;
			pp2[2] = r1.z;

			lambda = findIntersect(p1,p2,pp1,pp2,intersectpt);
			if (lambda >=0 && lambda <=1)
			{
				if (seg.m_tip[1].active == 1)
				{
					seg.m_tip[1].rt.x = intersectpt[0];
					seg.m_tip[1].rt.y = intersectpt[1];
					seg.m_tip[1].rt.z = intersectpt[2];
					seg.m_tip[1].active = 0;
				}
				else if (seg.m_tip[0].active == -1)
				{
					seg.m_tip[0].rt.x = intersectpt[0];
					seg.m_tip[0].rt.y = intersectpt[1];
					seg.m_tip[0].rt.z = intersectpt[2];
					seg.m_tip[0].active = 0;
				}
				cout << "3D intersection" << endl;
				++m_num_anastom;
			}
		}
	}


	return;
}


///////////////////////////////////////////////////////////////////////
// findIntersect
///////////////////////////////////////////////////////////////////////

// variables:
// a : start x,y point on segment1
// b : end x,y point on segment1
// c : start x,y point on segment2
// d : end x,y point on segment2
// t1 : displacement vector for segment1 (b-a)
// t2 : displacement vector for segment2 (d-c)
// lambda : distance to move along t1 from its start point
// mu : distance to move along t2 from its start point

// General case for the minimum distance between two lines
// d^2 = [BCD + B^2E + C^2F + A(D^2-4EF)]/[C^2 - 4AE]
// where
// A = t1.t1
// B = 2(a.t1 - t1.c)
// C = 2(t1.t2)
// D = 2(t2.c - t2.a)
// E = t2.t2
// F = a.a + c.c
// This will have a unique solution if
// | 2A -C |
// | -C	2E | is non-zero
// closest point on t1 is : a + lambda*t1
// where
// lambda = (C*mu - B)/(2A)
// and
// mu = (2*A*D + B*C)/(C^2 - 4*A*E) 


double Culture::findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3])
{
	double t1[3], t2[3];
	double A, B, C, D, E, F, min_dist;
	double lambda, mu;
	
	t1[0] = b[0] - a[0];
	t1[1] = b[1] - a[1];
	t1[2] = b[2] - a[2];

	t2[0] = d[0] - c[0];
	t2[1] = d[1] - c[1];
	t2[2] = d[2] - c[2];

	A = vec_dot(t1,t1);
	B = 2*(vec_dot(a,t1) - vec_dot(t1,c));
	C = 2*(vec_dot(t1,t2));
	D = 2*(vec_dot(t2,c) - vec_dot(t2,a));
	E = vec_dot(t2,t2);
	F = vec_dot(a,a)+ vec_dot(c,c);

	if (4*A*E-C*C != 0)
	{
		min_dist = (B*C*D + B*B*E + C*C*F + A*(D*D-4*E*F))/(C*C - 4*A*E);
		if (min_dist <= 7)
		{
			mu = (2*A*D + B*C)/(C*C - 4*A*E);
			lambda = (C*mu - B)/(2*A);
			intersectpt[0] = a[0] + lambda*t1[0];
			intersectpt[1] = a[1] + lambda*t1[1];
			intersectpt[2] = a[2] + lambda*t1[2];
		}
		else
		{
			return 1000;
		}
	}
	else
	{
		return 1000;
	}

	return lambda;
}


///////////////////////////////////////////////////////////////////////
// intersectPlane
///////////////////////////////////////////////////////////////////////

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


double N0[3], N1[3], N2[3], N3[3], N4[3], N5[3]; //face normals
double P0[3], P1[3], P2[3], P3[3], P4[3], P5[3]; //point on faces


bool Culture::intersectPlane(Segment &seg, int n, double intersectpt[3])
{
	Grid& grid = m_angio.grid;

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
	P1[0] = grid.xrange[1];
	P1[1] = 0;
	P1[2] = 0;

	//back
	N2[0] = 0;
	N2[1] = 1;
	N2[2] = 0;
	P2[0] = 0;
	P2[1] = grid.yrange[1];
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
	P4[2] = grid.zrange[1];

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


	if (seg.m_tip[1].active == 1)
	{
		V[0] = seg.m_tip[1].rt.x - seg.m_tip[0].rt.x;
		V[1] = seg.m_tip[1].rt.y - seg.m_tip[0].rt.y;
		V[2] = seg.m_tip[1].rt.z - seg.m_tip[0].rt.z;
		LP[0] = seg.m_tip[0].rt.x;
		LP[1] = seg.m_tip[0].rt.y;
		LP[2] = seg.m_tip[0].rt.z;
	}
	else
	{
		V[0] = seg.m_tip[0].rt.x - seg.m_tip[1].rt.x;
		V[1] = seg.m_tip[0].rt.y - seg.m_tip[1].rt.y;
		V[2] = seg.m_tip[0].rt.z - seg.m_tip[1].rt.z;
		LP[0] = seg.m_tip[1].rt.x;
		LP[1] = seg.m_tip[1].rt.y;
		LP[2] = seg.m_tip[1].rt.z;
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
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
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
				if (intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
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
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
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
				if (intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
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
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1])
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
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1])
					return true;
			}
			break;
	}
	return false;
}








///////////////////////////////////////////////////////////////////////
// PeriodicBC
///////////////////////////////////////////////////////////////////////

//table of faces opposite to 0,1,..6
double oppface[6];

Segment Culture::PeriodicBC(Segment &seg)
{
	Grid& grid = m_angio.grid;

	oppface[0] = grid.yrange[1];
	oppface[1] = grid.xrange[0];
	oppface[2] = grid.yrange[0];
	oppface[3] = grid.xrange[1];
	oppface[4] = grid.zrange[0];
	oppface[5] = grid.zrange[1];
	double unit_vec[3] = {0};
	double length = 0.0;
	double rem_length = 0.0;
	int n = 0;
	double intersectpt[3] = {0};
	
	if (seg.m_tip[1].active == 1)
	{
		unit_vec[0] = (seg.m_tip[1].rt.x-seg.m_tip[0].rt.x);
		unit_vec[1] = (seg.m_tip[1].rt.y-seg.m_tip[0].rt.y);
		unit_vec[2] = (seg.m_tip[1].rt.z-seg.m_tip[0].rt.z);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		
		for (n=0;n<6;++n)
		{
			if (intersectPlane(seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.label << " Crossed Face " << n << endl;
				seg.m_tip[1].rt.x = intersectpt[0];
				seg.m_tip[1].rt.y = intersectpt[1];
				seg.m_tip[1].rt.z = intersectpt[2];
				if (seg.length < 0)
				{
					seg.length = -(seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
				}
				else
				{
					seg.length = (seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
				}

				seg.m_tip[1].active = 0;
				seg.BCdead = 1;
				
				if (n > 3)
				    m_num_zdead += 1;
				return seg;
				
				m_angio.cult.m_frag.push_front (seg);
				Segment seg2;
				m_num_vessel = m_num_vessel + 1;
				
				switch (n)
				{
				case 0:
					seg2.m_tip[0].rt.x = seg.m_tip[1].rt.x;
					seg2.m_tip[0].rt.y = oppface[n];
					seg2.m_tip[0].rt.z = seg.m_tip[1].rt.z;
					break;
					
				case 1:
					seg2.m_tip[0].rt.x = oppface[n];
					seg2.m_tip[0].rt.y = seg.m_tip[1].rt.y;
					seg2.m_tip[0].rt.z = seg.m_tip[1].rt.z;
					break;
				case 2:
					seg2.m_tip[0].rt.x = seg.m_tip[1].rt.x;
					seg2.m_tip[0].rt.y = oppface[n];
					seg2.m_tip[0].rt.z = seg.m_tip[1].rt.z;
					break;
				case 3:
					seg2.m_tip[0].rt.x = oppface[n];
					seg2.m_tip[0].rt.y = seg.m_tip[1].rt.y;
					seg2.m_tip[0].rt.z = seg.m_tip[1].rt.z;
					break;
				case 4:
					seg2.m_tip[0].rt.x = seg.m_tip[1].rt.x;
					seg2.m_tip[0].rt.y = seg.m_tip[1].rt.y;
					seg2.m_tip[0].rt.z = oppface[n];
					break;
				case 5:
					seg2.m_tip[0].rt.x = seg.m_tip[1].rt.x;
					seg2.m_tip[0].rt.y = seg.m_tip[1].rt.y;
					seg2.m_tip[0].rt.z = oppface[n];
					break;
				}

				if (seg.length > 0) 
				{
					rem_length = length - (seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
					seg2.m_tip[1].rt.x = seg2.m_tip[0].rt.x + rem_length*unit_vec[0];
					seg2.m_tip[1].rt.y = seg2.m_tip[0].rt.y + rem_length*unit_vec[1];
					seg2.m_tip[1].rt.z = seg2.m_tip[0].rt.z + rem_length*unit_vec[2];
				}
				else 
				{
					rem_length = -length + (seg.m_tip[1].rt- seg.m_tip[0].rt).norm();
					seg2.m_tip[1].rt.x = seg2.m_tip[0].rt.x - rem_length*unit_vec[0];
					seg2.m_tip[1].rt.y = seg2.m_tip[0].rt.y - rem_length*unit_vec[1];
					seg2.m_tip[1].rt.z = seg2.m_tip[0].rt.z - rem_length*unit_vec[2];
				}
				
				seg2.m_tip[1].active = 1;
				seg2.m_tip[0].active = 0;
				seg2.label = seg.label;
				seg2.vessel = m_num_vessel;
				seg2.m_sprout = seg.m_sprout;
				seg2.length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				seg2.Recent_branch = seg.Recent_branch;
				seg2.TofBirth = seg.TofBirth;
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2);
				else
					return seg2;
			}
			
		}
	}
	
	else
	{
		unit_vec[0] = (seg.m_tip[0].rt.x-seg.m_tip[1].rt.x);
		unit_vec[1] = (seg.m_tip[0].rt.y-seg.m_tip[1].rt.y);
		unit_vec[2] = (seg.m_tip[0].rt.z-seg.m_tip[1].rt.z);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		for (n=0;n<6;++n)
		{
			if (intersectPlane(seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.label << " Crossed Face " << n << endl;
				seg.m_tip[0].rt.x = intersectpt[0];
				seg.m_tip[0].rt.y = intersectpt[1];
				seg.m_tip[0].rt.z = intersectpt[2];
				if (seg.length < 0)
				{
					seg.length = -(seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
				}
				else
				{
					seg.length = (seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
				}
				
				seg.m_tip[0].active = 0;
                
				if (n > 3)
				    m_num_zdead += 1;
				return seg;
				
				m_angio.cult.m_frag.push_front (seg);
				Segment seg2;
				m_num_vessel = m_num_vessel + 1;
				
				switch (n)
				{
				case 0:
					seg2.m_tip[1].rt.x = seg.m_tip[0].rt.x;
					seg2.m_tip[1].rt.y = oppface[n];
					seg2.m_tip[1].rt.z = seg.m_tip[0].rt.z;
					break;
					
				case 1:
					seg2.m_tip[1].rt.x = oppface[n];
					seg2.m_tip[1].rt.y = seg.m_tip[0].rt.y;
					seg2.m_tip[1].rt.z = seg.m_tip[0].rt.z;
					break;
				case 2:
					seg2.m_tip[1].rt.x = seg.m_tip[0].rt.x;
					seg2.m_tip[1].rt.y = oppface[n];
					seg2.m_tip[1].rt.z = seg.m_tip[0].rt.z;
					break;
				case 3:
					seg2.m_tip[1].rt.x = oppface[n];
					seg2.m_tip[1].rt.y = seg.m_tip[0].rt.y;
					seg2.m_tip[1].rt.z = seg.m_tip[0].rt.z;
					break;
				case 4:
					seg2.m_tip[1].rt.x = seg.m_tip[0].rt.x;
					seg2.m_tip[1].rt.y = seg.m_tip[0].rt.y;
					seg2.m_tip[1].rt.z = oppface[n];
					break;
				case 5:
					seg2.m_tip[1].rt.x = seg.m_tip[0].rt.x;
					seg2.m_tip[1].rt.y = seg.m_tip[0].rt.y;
					seg2.m_tip[1].rt.z = oppface[n];
					break;
				}

				if (seg.length < 0) 
				{
					rem_length = -length + (seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
					seg2.m_tip[0].rt.x = seg2.m_tip[1].rt.x - rem_length*unit_vec[0];
					seg2.m_tip[0].rt.y = seg2.m_tip[1].rt.y - rem_length*unit_vec[1];
					seg2.m_tip[0].rt.z = seg2.m_tip[1].rt.z - rem_length*unit_vec[2];
					
				}
				else 
				{
					rem_length = length - (seg.m_tip[1].rt - seg.m_tip[0].rt).norm();
					seg2.m_tip[0].rt.x = seg2.m_tip[1].rt.x + rem_length*unit_vec[0];
					seg2.m_tip[0].rt.y = seg2.m_tip[1].rt.y + rem_length*unit_vec[1];
					seg2.m_tip[0].rt.z = seg2.m_tip[1].rt.z + rem_length*unit_vec[2];
				}
				
				seg2.m_tip[1].active = 0;
				seg2.m_tip[0].active = -1;
				seg2.label = seg.label;
				seg2.vessel = m_num_vessel;
				seg2.m_sprout = seg.m_sprout;
				seg2.length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				seg2.Recent_branch = seg.Recent_branch;
				seg2.TofBirth = seg.TofBirth;
				
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2);
				else
					return seg2;
			}
		}
	}
	return seg;
}
