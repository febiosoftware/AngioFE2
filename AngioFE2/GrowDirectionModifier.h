#pragma once

#include "StdAfx.h"
#include "Segment.h"
class Culture;

//the interface for all operations that modify the growth direction
//some examples of this include deflection due to gradient and change in direction for anastamosis
//in the future this can be used to change the direction based on vegf concentrations
class GrowDirectionModifier
{
public:
	GrowDirectionModifier(Culture * c, int p) : culture(c), priority(p){}
	virtual ~GrowDirectionModifier(){}
	virtual vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) = 0;
	//used to sort these by priority
	int Priority() const { return priority; }
protected:
	Culture * culture;
	int priority;

};
//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction
class DefaultGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DefaultGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//this class changes the grow direction if the segment is a branch
class BranchGrowDirectionModifier : public GrowDirectionModifier
{
public:
	BranchGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};

//the class modifies the grow dierction if the gradeint is above a given threshold
class GradientGrowDirectionModifier : public GrowDirectionModifier
{
public:
	GradientGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};
//modifies the direction a segment grows based on its proximity to other segments
class AnastamosisGrowDirectionModifier : public GrowDirectionModifier
{
public:
	AnastamosisGrowDirectionModifier(Culture * c, int p) : GrowDirectionModifier(c, p){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, bool branch) override;
};