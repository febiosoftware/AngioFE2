#pragma once
#include <vector>
#include <FECore/vec3d.h>

class FEAngioNodeData
{
public:
	FEAngioNodeData();
	double m_ecm_den;		// current ecm density
	double m_ecm_den0;		// initial ecm density

	vec3d m_collfib;		// current collagen fiber direction
	vec3d m_collfib0;		// initial collagen fiber direction
	double vegf_conc;		// vegf concentration
	double m_da;			// ansiotropy value
	double alpha;

	//are id and n tag needed?
	//ntag is used in calculating the ecm density
	//in some way conpensates for multiple runs over the same node
	int m_ntag;
};


class FEAngioElementData
{
public:
	FEAngioElementData();
	double alpha;
	vec3d fiber_orient;
	std::vector<int> surfacesIndices;
};

//-----------------------------------------------------------------------------
inline double max3(double a, double b, double c)
{
	if ((a >= b)&&(a >= c)) return a;
	if ((b >= a)&&(b >= c)) return b;
	return c;
}