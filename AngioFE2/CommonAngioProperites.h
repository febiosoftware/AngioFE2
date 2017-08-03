#pragma once
#include <FECore/FEMaterial.h>
#include "FEProbabilityDistribution.h"
#include <FEBioMech/FESolidMaterial.h>
#include "FiberManager.h"
#include "FragmentBranching.h"
#include "BC.h"
#include "FragmentSeeder.h"
#include "GrowDirectionModifier.h"

class CommonAngioProperties :public FEMaterial
{
public:
	CommonAngioProperties(FEModel * pfem);
	~CommonAngioProperties(){}

	FEPropertyT<GrowDirectionModifiers> gdms;
	FEPropertyT<FragmentSeeder> fseeder;
	FEPropertyT<BC> bc;
	FEPropertyT<FragmentBranching> fbrancher;
	FEPropertyT<FESolidMaterial> vessel_material;
	FEPropertyT<FEProbabilityDistribution> length_to_branch;
	FEPropertyT<FiberInitializer> fiber_initializer;

	void InitializeFibers(FiberManager * man);
	void UpdateGDMs();
	void SetCulture(Culture * culture);
};
