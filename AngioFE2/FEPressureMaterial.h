#pragma once
#include <FEBioMech/FEElasticMaterial.h>

//-----------------------------------------------------------------------------
class FEPressureMaterial : public FEElasticMaterial
{
public:
	explicit FEPressureMaterial(FEModel* pfem) : FEElasticMaterial(pfem) { m_p = 0; }

	mat3ds Stress(FEMaterialPoint& mp) override;

	tens4ds Tangent(FEMaterialPoint& mp) override;

	double	m_p;	// pressure

	DECLARE_PARAMETER_LIST();
};