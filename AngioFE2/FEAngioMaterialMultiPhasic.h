#pragma once
#include "CommonAngioProperites.h"
#include "FEAngioMaterialBase.h"
#include <FEBioMix/FEMultiphasicStandard.h>

//#define N_MF

#ifndef N_MF
//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterialMultiPhasic : public FEMultiphasic, public FEAngioMaterialBase
{
public:


	explicit FEAngioMaterialMultiPhasic(FEModel* pfem);
	virtual ~FEAngioMaterialMultiPhasic();

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

	FEMaterial * GetMatrixMaterial() override { return matrix_material; }

	CommonAngioProperties * GetCommonAngioProperties() override { return common_properties; };

	int GetID_ang() const override { return FEMultiphasic::GetID(); };

	FEMaterial * GetMaterial()override { return dynamic_cast<FEMaterial*>(this); }

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

	//double StrainEnergyDensity(FEMaterialPoint& mp) override;

	// we use this to define a sprout in the material section of the input file
	void SetParameter(FEParam& p) override;

	void SetupSurface() override;

	FEElasticMaterial*			GetSolid() override { return matrix_material->GetSolid(); }
	FEHydraulicPermeability*	GetPermeability() override { return matrix_material->GetPermeability(); }
	FEOsmoticCoefficient*		GetOsmoticCoefficient() override { return  matrix_material->GetOsmoticCoefficient(); }
	FESolventSupply*			GetSolventSupply() override { return  matrix_material->GetSolventSupply(); }
	FESolute*					GetSolute(int i) override { return  matrix_material->GetSolute(i); }
	FESolidBoundMolecule*		GetSBM(int i) override { return  matrix_material->GetSBM(i); }
	FEChemicalReaction*			GetReaction(int i) override { return matrix_material->GetReaction(i); }


	int Solutes() override { return  matrix_material->Solutes(); }
	int SBMs() override { return matrix_material->SBMs(); }
	int Reactions() override { return matrix_material->Reactions(); }
private:
	DECLARE_PARAMETER_LIST();

	FEPropertyT<FEMultiphasic> matrix_material;
public:
	void ApplySym() override;
	FEPropertyT<CommonAngioProperties> common_properties;
	
};
#endif