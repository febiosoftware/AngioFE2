#pragma once

#include "StdAfx.h"
#include "Segment.h"
#include "FECore/FEMaterial.h"
#include <FECore/DataRecord.h>
#include <FECore/Archive.h>
#include "FEBioPlot/FEBioPlotFile2.h"
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
	//will be called once before growth per FE timestep
	virtual void Update() {}
	virtual vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) = 0;
	//used to sort these by priority

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);

protected:
	Culture * culture;
};

//
class GDMArchive : public Archive
{
public:
	GDMArchive(){}
	void reset();

	void WriteData(int nid, std::vector<float>& data) override;
	//domain,element is 0 indexed these work with per element data
	mat3ds GetDataMat3ds(int domain, int element_index);
	mat3dd GetDataMat3dd(int domain, int element_index);
	mat3d  GetDataMat3d(int domain, int element_index);
	float  GetDataFloat(int domain, int element_index);
	vec3d  GetDataVec3d(int domain, int element_index);
	//consider adding diagonal matrix and 4d tensor
	
private:
	std::vector<std::vector<float>> fpdata;
};

//a material which has a collection of grow direction modifiers
class GrowDirectionModifiers : public FEMaterial
{
public:
	GrowDirectionModifiers(FEModel * model);
	virtual ~GrowDirectionModifiers(){}

	vec3d ApplyModifiers(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length);

	void SetCulture(Culture * c);

	void Update();

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

//this class changes the segment length based on the average segment length load curve, cannot be the inital segment_length modifier
class DataStoreLengthDoubleGrowDirectionModifier : public GrowDirectionModifier
{
public:
	bool Init() override;
	DataStoreLengthDoubleGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;

private:
	int record_index;
	int field = 0;
	char field_name[DataRecord::MAX_STRING];
	DECLARE_PARAMETER_LIST();
};

//this class changes the segment length based on the average segment length load curve, cannot be the inital segment_length modifier
class PlotFile2DoubleGrowDirectionModifier : public GrowDirectionModifier
{
public:
	bool Init() override;
	void Update() override;
	PlotFile2DoubleGrowDirectionModifier(FEModel * model): GrowDirectionModifier(model){}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;

private:
	list<FEBioPlotFile2::DICTIONARY_ITEM>::const_iterator record_index;
	int field = 0;
	char field_name[DataRecord::MAX_STRING];
	//! means not currenlty supported
	int datatype = 0;//0 Float, 1 vec3d,2 Mat3ds, 3 Mat3dd!,4 tens4d!, 5 Mat3d
	GDMArchive archive;
	DECLARE_PARAMETER_LIST();
};


//pure virtual base class for Generic Growth Parameters
//all matrix indices will be 1 indexed
//conventions vectors will be read into col1 of the matrix
//doubles will be read into position 11 of the matrix
class GGP : public FEMaterial
{
public:
	GGP(FEModel * model) : FEMaterial(model) { AddProperty(&child, "child"); child.m_brequired = false; }
	virtual ~GGP() {}
	//will be called once before growth per FE timestep
	virtual void Update() { if (child) { child->Update(); } };
	virtual mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) { if (child) { return child->Operation(in, mat, tip); } return in; };

	//must be called before anything else is done but construction
	void SetCulture(Culture * cp) { culture = cp; }

protected:
	Culture * culture = nullptr;
	FEPropertyT<GGP> child;
};


class Plot2GGP : public GGP
{
public:
	Plot2GGP(FEModel * model) : GGP(model) {}
	virtual ~Plot2GGP() {}
	//will be called once before growth per FE timestep
	bool Init() override;

	void Update() override;
	mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	list<FEBioPlotFile2::DICTIONARY_ITEM>::const_iterator record_index;
	char field_name[DataRecord::MAX_STRING];
	//datatype is autmatically deteced 
	GDMArchive archive;
	DECLARE_PARAMETER_LIST();
};

class MatrixConverterGGP: public GGP
{
public:
	MatrixConverterGGP(FEModel * model);
	virtual ~MatrixConverterGGP() {}
	//will be called once before growth per FE timestep
	bool Init()override { return true; }

	void Update() override { GGP::Update(); }
	mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	double m[9][9];//the matrix that the unrolled matrix will be multiplied by
	DECLARE_PARAMETER_LIST();
};
//the nested tree will be executed before the child tree
class ForkedGGP : public GGP
{
public:
	//disallow empy forks
	ForkedGGP(FEModel * model) : GGP(model) { AddProperty(&child, "nest");}
	~ForkedGGP() {}
	//will be called once before growth per FE timestep
	void Update() override {};
	mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) override;


protected:
	Culture * culture = nullptr;
	FEPropertyT<GGP> nest;
};

//begin difficult to compute functions with the above classes

//computes the eigen values and stores them in the first collumn of the matrix all other values are 0
//must operate on a symmetric matrix if the in matrix is not symetric the results of this function is not defined
class EigenValuesGGP: public GGP
{
public:
	EigenValuesGGP(FEModel * model) : GGP(model) {  }
	virtual ~EigenValuesGGP() {}

	mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//computes the eigen vectors and stores them in the collums of the matrix
//must operate on a symmetric matrix if the in matrix is not symetric the results of this function is not defined
class EigenVectorsGGP : public GGP
{
public:
	EigenVectorsGGP(FEModel * model) : GGP(model) {  }
	virtual ~EigenVectorsGGP() {}

	mat3d Operation(mat3d in, FEAngioMaterial* mat, Segment::TIP& tip) override;
};