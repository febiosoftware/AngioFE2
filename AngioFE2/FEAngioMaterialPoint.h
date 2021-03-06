#pragma once
#include "StdAfx.h"
#include <FECore/FEMaterialPoint.h>
#include "Segment.h"

//-----------------------------------------------------------------------------
// A new material point class is defined to store the elastic parameters for 
// each integration point.
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	explicit FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt);

	//! The init function is used to intialize data
	void Init() override;

	//needed to properly evaluate viscoelastic materials
	void Update(const FETimeInfo& timeInfo) override;

	//! copy material point data (for running restarts) todo Is this still used?
	FEMaterialPoint* Copy() override;

	//! copy material point data (for running restarts) todo Is this still used?
	void Serialize(DumpStream& dmp) override;

	// These are the material parameters
	double		m_D;		// collagen density (?)
	double		m_DA;		// degree of anisotropy (?)
							//TODO: consider adding a custom weigth parameter per element 
	//the temporary to store angio stress
	mat3ds m_as;

	double ref_ecm_density;

	GridPoint	m_pt;	// grid point location of this material point

	double vessel_weight;
	double matrix_weight;
	FEMaterialPoint* vessPt;
	FEMaterialPoint* matPt;

	DECLARE_PARAMETER_LIST();

public:
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};