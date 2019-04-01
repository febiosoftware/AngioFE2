///////////////////////////////////////////////////////////////////////
// angio3d.h
///////////////////////////////////////////////////////////////////////


#pragma once

#include <math.h>

//TODO: consider replacing/removing most of this file

const double pi = 3.14159265;                                   // Global constant pi


///////////////////////// LINEAR ALGEBRA FUNCTIONS /////////////////////////

inline double vec_dot(double a[3], double b[3])
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double vec_norm(double a[3])
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

inline int sign(double n) { return (n>=0.0? 1 : -1); }


///////////////////////////// OTHER FUNCTIONS /////////////////////////////

inline double distance(double x0, double xf, double y0, double yf, double z0, double zf)
{
    return sqrt(pow(x0-xf,2.0)+pow(y0-yf,2.0)+pow(z0-zf,2.0));
}

// Insert function that constricts the angle into the domain +/- pi
inline double ang_dom(double angle)
{
    if (angle != 0.0)
    {
        while(1)
        {
            if (angle == pi)
                return angle;
            
            if (angle == -pi)
                return angle;
            
            if (angle > pi)
                angle = angle - 2*pi;
    
            if (angle < -pi)
                angle = angle + 2*pi;
        
            if ((angle <= pi) && (angle >= -pi))
                return angle;
        }
    }
    
    return angle;
}

// generates a random number between 0, 1
inline double frand() { return (double) rand() / RAND_MAX; }

// generates a random vector in the unit cube [-1, +1]
//TODO: This doesn't generate a truly random vector since there will be bias towards corners.
//when generating a direction with this function the vectors will be ~30% different from a uniform distribution of
//directions
inline vec3d vrand() { return vec3d(2.0*(frand()-0.5), 2.0*(frand()-0.5), 2*(frand() - 0.5)); }

//see https://www.opengl.org/sdk/docs/man/html/reflect.xhtml
inline vec3d reflect(vec3d & I, vec3d & N)
{
	double temp = (2.0 * (N * I));
	return I - (N * temp);
}
//see https://www.opengl.org/sdk/docs/man/html/mix.xhtml
inline vec3d mix(vec3d & x, vec3d & y, double a)
{
	return x*(1 - a) + y*a;
}
//consider making this a template or using GLM
inline double mix(double x, double y, double a)
{
	return x*(1 - a) + y*a;
}
//new mix method 
//x is per, y is col_dir
// method of rotation about an axis: https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
inline vec3d mix3d(vec3d & x, vec3d & y, double a)
{
	// determine the angle between the two vectors then scale it
	double phi = a * acos(x * y);
	// determine the normal vector to the plane spanned by x and y
	vec3d normal_vec = x ^ y;
	normal_vec.unit();
	if (normal_vec.z < 0) phi = -phi;
	mat3d rot_mat;
	auto nx = normal_vec.x, ny = normal_vec.y, nz = normal_vec.z;
	double cp = cos(phi), sp = sin(phi);
	// assemble the rotation matrix about the normal vector within the plane spanned by x and y
	rot_mat[0][0] = cp + pow(nx,2) * (1 - cp);	rot_mat[0][1] = nx * ny*(1 - cp) - nz * sp;	rot_mat[0][2] = nx * nz*(1 - cp) + ny * sp;
	rot_mat[1][0] = ny * nx*(1 - cp) + nz * sp;	rot_mat[1][1] = cp + pow(ny, 2)*(1 - cp);	rot_mat[1][2] = ny * nz*(1 - cp) - nx * sp;
	rot_mat[2][0] = nz * nx*(1 - cp) - ny * sp;	rot_mat[2][1] = nz * ny*(1 - cp) + nx * sp;	rot_mat[2][2] = cp + pow(nz, 2)*(1 - cp);
	return rot_mat * x;
}

//binary search used in volume/area seeders
static size_t findElement(double val, int lo, int high, double * begin, double * end)
{
	int mid = lo + (high - lo) / 2;
	if (val < begin[mid])
	{
		return findElement(val, lo, mid - 1, begin, end);
	}
	else if (val > end[mid])
	{
		return findElement(val, mid + 1, high, begin, end);
	}
	else
	{
		return mid;
	}
}