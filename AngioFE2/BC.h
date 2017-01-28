#pragma once
#include "StdAfx.h"
#include <FECore/vec3d.h>
#include "FECore/FEMaterial.h"


//-----------------------------------------------------------------------------
class FEAngio;
class Elem;
class Segment;
class Node;
class FESolidElement;
class Culture;
class FEAngioMaterial;


//multAngioBC is boundary condition used when a segment would cross from one angio material to another
class MBC : public FEMaterial
{
public:
	virtual ~MBC(){}
	MBC(FEModel * model);
	//returns whether or not this class is going to hanel the boundary between the two angio materials
	virtual bool acceptBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1);
	virtual void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) = 0;

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);
protected:
	Culture * culture;
};
class SameMBC : public MBC
{
public:
	SameMBC(FEModel * model);
	~SameMBC(){}
	//returns whether or not this class is going to hanel the boundary between the two angio materials
	void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) override { assert(false); }
};


//-----------------------------------------------------------------------------
// The BC class is used by the CULTURE class to handle boundary 
// conditions within the model.
class BC : public FEMaterial
{
public:
	BC(FEModel * model);

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);

	virtual ~BC();
	
	// checks if a new segment has cross the boundary
	//tip 0 is assumed to be in bounds
	//returns the address of the added segment
	void CheckBC(Segment &seg);

	//checks if a segment would leave it's current meterial
	//returns true if the material changes
	//used in anastomosis
	bool ChangeOfMaterial(Segment & seg) const;

	const double epsilon = 0.0001;//any elements with any dimension smaller than this may not be properly handled by collision detection
	//consider looking at nextafter and nexttoward to replace this. the difficulty is they wil be in terms of the whole vector

protected:
	virtual void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) = 0;
	Culture * culture= nullptr;
	FEPropertyT<MBC> mbc;
private:
	BC & operator=(const BC&);
};

class BouncyBC: public BC
{
public:
	BouncyBC(FEModel * model);
	virtual ~BouncyBC(){}
protected:
	void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	BouncyBC & operator=(const BouncyBC&);
};
class StopBC: public BC
{
public:
	StopBC(FEModel * model);
	virtual ~StopBC(){}
protected:
	void HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se) override;
private:
	StopBC & operator=(const StopBC&);
};
class PassThroughMBC : public MBC
{
public:
	PassThroughMBC(FEModel * model);
	~PassThroughMBC(){}
	void handleBoundary(FEAngioMaterial * mat0, FEAngioMaterial * mat1, Segment & seg) override;
protected:
};