#pragma once

#include "StdAfx.h"
#include "Segment.h"
#include "FECore/FEMaterial.h"
#include "FEProbabilityDistribution.h"

class FEAngioMaterialBase;
class FiberManager;

//the base class for classes that initialize the fibers
class FiberInitializer : public FEMaterial
{
public:
	explicit FiberInitializer(FEModel * model):FEMaterial(model){}
	virtual ~FiberInitializer(){}
	virtual void InitializeFibers(FiberManager * fman) = 0;

	void nodeToInt(FiberManager * fman);
	void GridPointOfIntPoint(FESolidElement * se, int ei, int intp, GridPoint & gp);
};
//does nothing the fiber direction will be unchanged
class NullFiberInitializer : public FiberInitializer
{
public:
	explicit NullFiberInitializer(FEModel * model) : FiberInitializer(model){}
	virtual ~NullFiberInitializer(){}
	void InitializeFibers(FiberManager * fman) override{}
};
//set the fibers to a random orientation
class RandomFiberInitializer : public FiberInitializer
{
public:
	explicit RandomFiberInitializer(FEModel * model) : FiberInitializer(model){}
	virtual ~RandomFiberInitializer(){}
	void InitializeFibers(FiberManager * fman) override;
};

//set the fibers to a random orientation without mangling the relationship between the material axes at the integration points
class RandomFiberInitializerNonMangling : public FiberInitializer
{
public:
	explicit RandomFiberInitializerNonMangling(FEModel * model) : FiberInitializer(model) {}
	virtual ~RandomFiberInitializerNonMangling() {}
	void InitializeFibers(FiberManager * fman) override;
};
//set the fibers to a random orientation on a per element basis
class RandomFiberInitializerPE : public FiberInitializer
{
public:
	explicit RandomFiberInitializerPE(FEModel * model) : FiberInitializer(model) {}
	virtual ~RandomFiberInitializerPE() {}
	void InitializeFibers(FiberManager * fman) override;
};
class ExplicitDistributionsFiberInitializer: public FiberInitializer
{
public:
	explicit ExplicitDistributionsFiberInitializer(FEModel * model);
	virtual ~ExplicitDistributionsFiberInitializer() {}
	void Setup();
	void InitializeFibers(FiberManager * fman) override;
private:
	FEPropertyT<FEProbabilityDistribution> alpha;
	FEPropertyT<FEProbabilityDistribution> beta;
	FEPropertyT<FEProbabilityDistribution> gamma;
};

class EllipsoidalFiberInitializer : public FiberInitializer
{
public:
	explicit EllipsoidalFiberInitializer(FEModel * model);
	virtual ~EllipsoidalFiberInitializer()
	{
		if (totalWeightsBegin)
			delete[] totalWeightsBegin;
		if (totalWeightsEnd)
			delete[] totalWeightsEnd;
		if (directions)
			delete[] directions;
	}
	void InitializeFibers(FiberManager * fman) override;
private:
	DECLARE_PARAMETER_LIST();
	double a=2,b=1,c=1;
	int theta_slice = 360;
	int phi_slice = 360;
	double * totalWeightsBegin = nullptr;
	double * totalWeightsEnd = nullptr;
	vec3d * directions = nullptr;
};


class FiberManager
{
public:
	explicit FiberManager(FEAngioMaterialBase * mat) : material(mat){  }
	virtual ~FiberManager(){}

	vec3d GetFiberDirection(GridPoint & pt, double& lambda);
	vec3d GetMinorAxisDirection1(GridPoint & pt, double  &lambda);
	vec3d GetMinorAxisDirection2(GridPoint & pt, double & lambda);
	vec3d GetFiberAtNode(int node, double & lambda);
	vec3d GetMinor1AtNode(int node, double  &lambda);
	vec3d GetMinor2AtNode(int node, double  &lambda);
	void Update();

private:
	FEAngioMaterialBase * material;
	std::vector<std::vector<double>> fiber_at_int_pts[4];
	std::vector<std::vector<double>> m1_at_int_pts[4];
	std::vector<std::vector<double>> m2_at_int_pts[4];
	std::vector<double> fibers_at_nodes[4];
	std::vector<double> minoraxis1_at_nodes[4];
	std::vector<double> minoraxis2_at_nodes[4];

	friend class FiberInitializer;
	friend class RandomFiberInitializer;
	friend class RandomFiberInitializerNonMangling;
	friend class RandomFiberInitializerPE;
	friend class ExplicitDistributionsFiberInitializer;
	friend class EllipsoidalFiberInitializer;
};

