#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include <FEBioMech/FEElasticFiberMaterial.h>
#include <FECore/FEDataArray.h>
#include <FECore/FESurface.h>
#include <FECore/FENormalProjection.h>
#include "FEAngio.h"
#include "Culture.h"
#include "FEProbabilityDistribution.h"
#include "KDTree/kdtree.h"
#include "FiberManager.h"
#include "FEAngioMaterialPoint.h"
#include <FEBioMix/FEBiphasic.h>
#include "ECMInitializer.h"
#include "FEAngioMaterialBase.h"



//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterial : public FEElasticFiberMaterial, public FEAngioMaterialBase
{
public:
	

	explicit FEAngioMaterial(FEModel* pfem);
	virtual ~FEAngioMaterial();

	friend class Fileout;
	friend class FEPlotMatrixStress;
	friend class FEPlotVesselStress;
	friend class FEPlotMatrixTangent;

	//begin functions from FEAngioMaterialBase

	// Calculate the active Angio stress
	mat3ds AngioStress(FEAngioMaterialPoint& mp) override;

	void FinalizeInit() override;

	void UpdateECM() override;

	void UpdateGDMs() override;

	bool InitECMDensity(FEAngio * angio)  override;

	void InitializeFibers() override;
	
	//begin functions from FEMaterial

	// material initialization
	bool Init() override;

	// Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp) override;


	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override;

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

	double StrainEnergyDensity(FEMaterialPoint& mp) override;

	// we use this to define a sprout in the material section of the input file
	void SetParameter(FEParam& p) override;

	double GetAnisotropy() const;
	
	void SetupSurface();
private:
	DECLARE_PARAMETER_LIST();


public:
	void ApplySym();
	FEPropertyT<GrowDirectionModifiers> gdms;
	FEPropertyT<FragmentSeeder> fseeder;
	FEPropertyT<BC> bc;
	FEPropertyT<FragmentBranching> fbrancher;
	FEPropertyT<FESolidMaterial> vessel_material;
	FEPropertyT<FEMaterial> matrix_material;
private:
	FEPropertyT<FEProbabilityDistribution> length_to_branch;
	FEPropertyT<FiberInitializer> fiber_initializer;
};

//-----------------------------------------------------------------------------
class FEPressureMaterial : public FEElasticMaterial
{
public:
	explicit FEPressureMaterial(FEModel* pfem) : FEElasticMaterial(pfem){ m_p = 0; }

	mat3ds Stress(FEMaterialPoint& mp) override;

	tens4ds Tangent(FEMaterialPoint& mp) override;

	double	m_p;	// pressure

	DECLARE_PARAMETER_LIST();
};