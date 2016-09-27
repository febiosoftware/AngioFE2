#pragma once
#include <FECore/vec3d.h>

//-----------------------------------------------------------------------------
class FEAngio;
class Elem;
class Segment;
class Node;
class FESolidElement;
class Culture;

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

	BC(FEAngio& angio);
	virtual ~BC();
	
	// checks if a new segment has cross the boundary
	//tip 0 is assumed to be in bounds
	void CheckBC(Segment &seg, Culture * culture);

	//checks if a segment would leave it's current meterial
	//returns true if the material changes
	//used in anastomosis
	bool ChangeOfMaterial(Segment & seg) const;

	const double epsilon = 0.0001;//any elements with any dimension smaller than this may not be properly handled by collision detection
	//consider looking at nextafter and nexttoward to replace this. the difficulty is they wil be in terms of the whole vector

protected:
	virtual void HandleBoundary(Segment & seg, Culture * culture, vec3d lastGoodPt, double * rs, FESolidElement * se) = 0;
	FEAngio&	m_angio;
private:
	BC & operator=(const BC&);
};

class BouncyBC: public BC
{
public:
	BouncyBC(FEAngio & angio): BC(angio){}
	virtual ~BouncyBC(){}
protected:
	void HandleBoundary(Segment & seg, Culture * culture, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	BouncyBC & operator=(const BouncyBC&);
};
class StopBC: public BC
{
public:
	StopBC(FEAngio & angio) : BC(angio){}
	virtual ~StopBC(){}
protected:
	void HandleBoundary(Segment & seg, Culture * culture, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	StopBC & operator=(const StopBC&);
};