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
	virtual void SetCulture(Culture * cp);

protected:
	Culture * culture = nullptr;
};

//an archive that can store a plot2 variable from the previous finite element iteration
//this data is stored in ram. This class is used in GGPPlot2 classes
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
	//per node versions of the functions 
	mat3ds GetDataMat3ds(int domain, int element_index, Segment::TIP& tip);
	mat3dd GetDataMat3dd(int domain, int element_index, Segment::TIP& tip);
	mat3d  GetDataMat3d(int domain, int element_index, Segment::TIP& tip);
	float  GetDataFloat(int domain, int element_index, Segment::TIP& tip);
	vec3d  GetDataVec3d(int domain, int element_index, Segment::TIP& tip);
	//consider adding diagonal matrix and 4d tensor
	
	//gradient version of elemnet functions
	vec3d  GetDataGradientFloat(int domain, int element_index, Segment::TIP& tip, int size, int offset);
	void GradientEnabled(bool st) { gradient_defined = st; }
private:
	std::vector<std::vector<float>> fpdata;
	std::vector<double> unrolled_data;
	bool gradient_defined = false;
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


//pure virtual base class for Generic Growth Parameters
//all matrix indices will be 1 indexed
//conventions vectors will be read into diagonal of the matrix all other values will be zeroed
//conventions doubles will be read into element 11 of the matrix all other values will be zeroed
//child should not be required for any operation and should be the last operator applied by a GGP
class GGP : public FEMaterial
{
public:
	GGP(FEModel * model) : FEMaterial(model) { AddProperty(&child, "child"); child.m_brequired = false; }
	virtual ~GGP() {}
	//will be called once before growth per FE timestep
	virtual void Update() { if (child) { child->Update(); } };
	virtual mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) { if (child) { return child->Operation(in, fin, mat, tip); } return in; };

	bool Init() override
	{
		if (child)
			return child->Init();
		return true;
	}
	//must be called before anything else is done but construction
	virtual void SetCulture(Culture * cp)
	{
		culture = cp;
		if(child)
		{
			child->SetCulture(cp);
		}
	}

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
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

protected:
	list<FEBioPlotFile2::DICTIONARY_ITEM>::const_iterator record_index;
	char field_name[DataRecord::MAX_STRING];
	//datatype is autmatically deteced 
	GDMArchive archive;
	DECLARE_PARAMETER_LIST();
};

class GradientPlot2GGP : public Plot2GGP
{
public:
	GradientPlot2GGP(FEModel * model) : Plot2GGP(model) { archive.GradientEnabled(true); }
	virtual ~GradientPlot2GGP() {}
	//will be called once before growth per FE timestep
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	int size=1;
	int offset = 0;
	DECLARE_PARAMETER_LIST();
};

class NodalDataGGP : public GGP
{
public:
	NodalDataGGP(FEModel * model) : GGP(model) { field_name[0] = 0; }
	virtual ~NodalDataGGP() {}

	void Update() override;
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

	void SetCulture(Culture * cp) override;
private:
	char field_name[DataRecord::MAX_STRING];
	int offset = 0;
	vector<double> data;
	DECLARE_PARAMETER_LIST();
};

class MatrixConverterGGP: public GGP
{
public:
	MatrixConverterGGP(FEModel * model);
	virtual ~MatrixConverterGGP() {}
	//will be called once before growth per FE timestep

	void Update() override { GGP::Update(); }
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	double m[9][9];//the matrix that the unrolled matrix will be multiplied by
	DECLARE_PARAMETER_LIST();
};

//the nested tree will be executed before the child tree
class ForkedGGP : public GGP
{
public:
	//disallow empy forks
	ForkedGGP(FEModel * model) : GGP(model) { AddProperty(&nest, "nest"); }
	~ForkedGGP() {}
	//will be called once before growth per FE timestep
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
	void Update() override { 
		GGP::Update();
		nest->Update();
	}
	void SetCulture(Culture * cp) override
	{
		nest->SetCulture(cp);
		GGP::SetCulture(cp);
	}
	bool Init() override
	{
		if (nest->Init())
			return GGP::Init();
		return false;
	}

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

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//computes the eigen vectors and stores them in the collums of the matrix
//must operate on a symmetric matrix if the in matrix is not symetric the results of this function is not defined
class EigenVectorsGGP : public GGP
{
public:
	EigenVectorsGGP(FEModel * model) : GGP(model) {  }
	virtual ~EigenVectorsGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//the computes the cross product and stores it in the diagonal of the returned matrix
class CrossGGP : public GGP
{
public:
	//disallow empy forks
	CrossGGP(FEModel * model) : GGP(model)
	{
		AddProperty(&v1, "v1");
		AddProperty(&v2, "v2");
	}
	~CrossGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
	void Update() override {
		GGP::Update();
		v1->Update();
		v2->Update();
	}
	void SetCulture(Culture * cp) override
	{
		v1->SetCulture(cp);
		v2->SetCulture(cp);
		GGP::SetCulture(cp);
	}
	bool Init() override
	{
		if (v1.Init() && v2.Init())
			return GGP::Init();
		return false;
	}

protected:
	Culture * culture = nullptr;
	FEPropertyT<GGP> v1;
	FEPropertyT<GGP> v2;
};

//if (threshold * (1,0,0)).x > 0 then statement will be executed and modify the result of child
class ThresholdGGP : public GGP
{
public:
	//disallow empy forks
	ThresholdGGP(FEModel * model) : GGP(model)
	{
		AddProperty(&statement, "statement");
		AddProperty(&threshold, "threshold");
	}
	~ThresholdGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
	void Update() override {
		GGP::Update();
		statement->Update();
		threshold->Update();
	}
	void SetCulture(Culture * cp) override
	{
		statement->SetCulture(cp);
		threshold->SetCulture(cp);
		GGP::SetCulture(cp);
	}
	bool Init() override
	{
		if (statement.Init() && threshold.Init())
			return GGP::Init();
		return false;
	}

protected:
	Culture * culture = nullptr;
	FEPropertyT<GGP> statement;
	FEPropertyT<GGP> threshold;
	vec3d vec = vec3d(1, 0, 0);
	double condition = 0.0;
	DECLARE_PARAMETER_LIST();
};

//applies arccos to all elements of the matrix
class ArcCosGGP : public GGP
{
public:
	ArcCosGGP(FEModel * model) : GGP(model) {  }
	virtual ~ArcCosGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//applies arccos to all elements of the matrix
class ArcSinGGP : public GGP
{
public:
	ArcSinGGP(FEModel * model) : GGP(model) {  }
	virtual ~ArcSinGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//applies cos to all elements of the matrix
class CosGGP : public GGP
{
public:
	CosGGP(FEModel * model) : GGP(model) {  }
	virtual ~CosGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//applies sin to all elements of the matrix
class SinGGP : public GGP
{
public:
	SinGGP(FEModel * model) : GGP(model) {  }
	virtual ~SinGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//calculates the inverse of in and passes this on
class MatrixInverseGGP : public GGP
{
public:
	MatrixInverseGGP(FEModel * model) : GGP(model) {  }
	virtual ~MatrixInverseGGP() {}

	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;
};

//GGP classes used in testing

class SetterGGP : public GGP
{
public:
	//the user must do all initialization themselves
	SetterGGP(FEModel * model): GGP(model)
	{
		invec = vec3d(1, 0, 0);
		m = mat3d(1,0,0,
			0,1,0,
			0,0,1);
	}
	virtual ~SetterGGP() {}
	//will be called once before growth per FE timestep

	void Update() override { GGP::Update(); }
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	mat3d m;//the matrix 
	vec3d invec;
	DECLARE_PARAMETER_LIST();
};

class AssertGGP : public GGP
{
public:
	//the user must do all initialization themselves
	AssertGGP(FEModel * model) : GGP(model) {
		m = mat3d(1, 0, 0,
			0, 1, 0,
			0, 0, 1);
	}
	virtual ~AssertGGP() {}
	//will be called once before growth per FE timestep

	void Update() override { GGP::Update(); }
	mat3d Operation(mat3d in, vec3d fin, FEAngioMaterial* mat, Segment::TIP& tip) override;

private:
	mat3d m;//the matrix
	double tolerance= 0.01;
	DECLARE_PARAMETER_LIST();
};

//begin bind points
//will ignore the previous direction and generate the direction a segmetn should grow based on collagen direction
class DefaultGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DefaultGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
	void Update() override;
	void SetCulture(Culture * cp) override;
private:
	FEPropertyT<GGP> collagen_direction;
	FEPropertyT<GGP> previous_direction;
	FEPropertyT<GGP> weight_interpolation;
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
	void Update() override;
	void SetCulture(Culture * cp) override;
private:
	FEPropertyT<GGP> collagen_direction;
	FEPropertyT<GGP> previous_direction;
};

//this sets the segment length to 1
class UnitLengthGrowDirectionModifier : public GrowDirectionModifier
{
public:
	UnitLengthGrowDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
	bool Init() override
	{
		if (length_modifier)
			return length_modifier->Init();
		return true;
	}

	void Update() override
	{
		if(length_modifier)
		{
			length_modifier->Update();
		}
	}
	void SetCulture(Culture * cp) override
	{
		if(length_modifier)
		{
			length_modifier->SetCulture(cp);
		}
	}

private:
	FEPropertyT<GGP> length_modifier;
};

//this class modifies the segment length by the density factor
class DensityScaleGrowDirectionModifier : public GrowDirectionModifier
{
public:
	DensityScaleGrowDirectionModifier(FEModel * model);
	bool Init() override
	{
		if(density_scale)
			return density_scale->Init();
		return true;
	}
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
	void Update() override;
	void SetCulture(Culture * cp) override;
private:
	FEPropertyT<GGP> density_scale;
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
	void Update() override;
	void SetCulture(Culture * cp) override;
private:
	double threshold = 0.01;//the threshold over which vessels will deflect on the gradient
	FEPropertyT<GGP> threshold_scale;
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

class EdgeDeflectorGrowDirectionModifier : public GrowDirectionModifier
{
public:
	EdgeDeflectorGrowDirectionModifier(FEModel * model);

	void SetCulture(Culture * cp) override;
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterial* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
	//used to sort these by priority
private:
	std::unordered_map<int,bool> edge_element;
};