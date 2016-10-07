#pragma once
#include <FECore/vec3d.h>
#include <cassert>

//-----------------------------------------------------------------------------
class FEAngio;
class Elem;
class Segment;
class Node;
class FESolidElement;
class Culture;
class FEAngioMaterial;


//multAngioBC is boundary condition used when a segment would cross from one angio material to another
class MBC
{
public:
	virtual ~MBC(){}
	MBC(Culture * c){ culture = c; }
	//returns whether or not this class is going to hanel the boundary between the two angio materials
	virtual bool acceptBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1) = 0;
	virtual void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) = 0;
protected:
	Culture * culture;
};
class SameMBC : public MBC
{
public:
	SameMBC(Culture * c): MBC(c){}
	~SameMBC(){}
	//returns whether or not this class is going to hanel the boundary between the two angio materials
	bool acceptBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1) override { return false; }
	void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) override { assert(false); }
};


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

	BC(FEAngio& angio, Culture * c);
	virtual ~BC();
	
	// checks if a new segment has cross the boundary
	//tip 0 is assumed to be in bounds
	void CheckBC(Segment &seg);

	//checks if a segment would leave it's current meterial
	//returns true if the material changes
	//used in anastomosis
	bool ChangeOfMaterial(Segment & seg) const;

	const double epsilon = 0.0001;//any elements with any dimension smaller than this may not be properly handled by collision detection
	//consider looking at nextafter and nexttoward to replace this. the difficulty is they wil be in terms of the whole vector

protected:
	virtual void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) = 0;
	FEAngio&	m_angio;
	Culture * culture= nullptr;
private:
	BC & operator=(const BC&);
	MBC * mbc= nullptr;
};

class BouncyBC: public BC
{
public:
	BouncyBC(FEAngio & angio, Culture * c): BC(angio, c){}
	virtual ~BouncyBC(){}
protected:
	void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	BouncyBC & operator=(const BouncyBC&);
};
class StopBC: public BC
{
public:
	StopBC(FEAngio & angio, Culture * c) : BC(angio, c){}
	virtual ~StopBC(){}
protected:
	void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	StopBC & operator=(const StopBC&);
};
class PassThroughMBC : public MBC
{
public:
	PassThroughMBC(Culture * c) : MBC(c){}
	~PassThroughMBC(){}
	bool acceptBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1) override { return true; }
	void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) override;
protected:
};