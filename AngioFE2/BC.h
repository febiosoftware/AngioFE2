#pragma once
#include <list>
using namespace std;
#include <FECore/vec3d.h>
#include "Grid.h"

//-----------------------------------------------------------------------------
class FEAngio;
class Elem;
class Segment;
class Node;

//-----------------------------------------------------------------------------
// The BC class is used by the CULTURE class to handle boundary 
// conditions within the model.
class BC
{
public:
	enum {
		STOP,		// the segment stops at the boundary
		BOUNCY		// bouncy wall
	};

public:
	BC(FEAngio& angio);
	~BC();
	
	// checks if a new segment has cross the boundary
	void CheckBC(Segment &seg);

private:
	// enforces a boundary condition
	void EnforceBC(Segment &seg, FACE_INTERSECTION& ic);

	// apply BC where vessel stops growing after hitting boundary
	void BCStop(Segment &seg, FACE_INTERSECTION& ic);

	// apply BC where the vessel bounces off the boundary
	void BCBouncy(Segment &seg, FACE_INTERSECTION& ic);

private:
	FEAngio&	m_angio;
};
