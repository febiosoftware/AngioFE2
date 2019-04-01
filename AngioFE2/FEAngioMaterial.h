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
#include "CommonAngioProperites.h"
#include <FEBioMix/FEMultiphasic.h>


//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterial : public FEElasticMaterial, public FEAngioMaterialBase
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

	void UpdateAngioStresses() override;

	bool InitECMDensity(FEAngio * angio)  override;

	void InitializeFibers() override;
	
	FEMaterial * GetMatrixMaterial() override { return matrix_material; }

	CommonAngioProperties * GetCommonAngioProperties() override { return common_properties; };

	int GetID_ang() override { 
		FECoreBase  * base = GetParent();
		if(base)
		{
			return base->GetID();
		}
		return FEElasticMaterial::GetID(); 
	};

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

	double StrainEnergyDensity(FEMaterialPoint& mp) override;

	// we use this to define a sprout in the material section of the input file
	void SetParameter(FEParam& p) override;
	
	void SetupSurface() override;
private:
	DECLARE_PARAMETER_LIST();

	FEPropertyT<FESolidMaterial> matrix_material;
	FEPropertyT<CommonAngioProperties> common_properties;
public:
	void ApplySym() override;
	
	
	
};