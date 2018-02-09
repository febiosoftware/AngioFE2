#pragma once
#include "StdAfx.h"
#include <FECore/vec3d.h>

class FEAngioNodeData
{
public:
	FEAngioNodeData();
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
	std::vector<int> surfacesIndices;
	int flags;//used for boolean flags bit 0 is segment crossed entire material
};

//-----------------------------------------------------------------------------
inline double max3(double a, double b, double c)
{
	if ((a >= b)&&(a >= c)) return a;
	if ((b >= a)&&(b >= c)) return b;
	return c;
}