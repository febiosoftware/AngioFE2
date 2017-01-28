#pragma once

#include "StdAfx.h"
#include "Segment.h"
#include "FECore/FEMaterial.h"
class Culture;

//the interface for all operations that modify the growth direction
//some examples of this include deflection due to gradient and change in direction for anastamosis
//in the future this can be used to change the direction based on vegf concentrations, also consider a modifier for the weight of the previous direction
class GrowDirectionModifier : public FEMaterial
{
public:
	GrowDirectionModifier(FEModel * model);
	virtual ~GrowDirectionModifier(){}
	virtual vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) = 0;
	//used to sort these by priority

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);

protected:
	Culture * culture;
};
//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction
class DefaultGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DefaultGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//this class changes the grow direction if the segment is a branch
class BranchGrowDirectionModifier : public GrowDirectionModifier
{
public:
	BranchGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//the class modifies the grow dierction if the gradeint is above a given threshold
class GradientGrowDirectionModifier : public GrowDirectionModifier
{
public:
	GradientGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;

private:
	double threshold = 0.01;//the threshold over which vessels will deflect on the gradient
	DECLARE_PARAMETER_LIST();
};
//modifies the direction a segment grows based on its proximity to other segments
class AnastamosisGrowDirectionModifier : public GrowDirectionModifier
{
public:
	AnastamosisGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;

private:
	double search_radius;
};