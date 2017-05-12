#pragma once

#include "StdAfx.h"
#include "Segment.h"
#include "FECore/FEMaterial.h"
class FEAngioMaterial;
class Culture;




//the interface for all operations that modify the growth direction
//some examples of this include deflection due to gradient and change in direction for anastamosis
//in the future this can be used to change the direction based on vegf concentrations, also consider a modifier for the weight of the previous direction
class GrowDirectionModifier : public FEMaterial
{
public:
	GrowDirectionModifier(FEModel * model);
	virtual ~GrowDirectionModifier(){}
	virtual vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) = 0;
	//used to sort these by priority

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);

protected:
	Culture * culture;
};

//a material which has a collection of grow direction modifiers
class GrowDirectionModifiers : public FEMaterial
{
public:
	GrowDirectionModifiers(FEModel * model);
	virtual ~GrowDirectionModifiers(){}

	vec3d ApplyModifiers(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length);

	void SetCulture(Culture * c);

private:
	FEVecPropertyT<GrowDirectionModifier> grow_direction_modifiers;
	Culture * culture;
};

//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction
class DefaultGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DefaultGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};

//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction and current stretch
class BaseFiberAwareGrowDirectionModifier : public GrowDirectionModifier
{
public:
	BaseFiberAwareGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};

//this class changes the grow direction if the segment is a branch
class BranchGrowDirectionModifier : public GrowDirectionModifier
{
public:
	BranchGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};

//this sets the segment length to 1
class UnitLengthGrowDirectionModifier : public GrowDirectionModifier
{
public:
	UnitLengthGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};

//this class modifies the segment length by the density factor
class DensityScaleGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DensityScaleGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};
//this class changes the segment length based on the average segment length load curve, cannot be the inital segment_length modifier
class SegmentLengthGrowDirectionModifier : public GrowDirectionModifier
{
public:
	SegmentLengthGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
};

//the class modifies the grow dierction if the gradeint is above a given threshold
class GradientGrowDirectionModifier : public GrowDirectionModifier
{
public:
	GradientGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;

private:
	double threshold = 0.01;//the threshold over which vessels will deflect on the gradient
	DECLARE_PARAMETER_LIST();
};
//modifies the direction a segment grows based on its proximity to other segments if a tip is within the radius specified the vessels direction will be set to grow towards that segment
class AnastamosisGrowDirectionModifier : public GrowDirectionModifier
{
public:
	AnastamosisGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;

private:
	double search_radius = 100.0;
	double search_multiplier = 1.0;
	DECLARE_PARAMETER_LIST();
};