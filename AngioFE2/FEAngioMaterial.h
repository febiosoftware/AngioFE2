#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include "FEAngio.h"

//-----------------------------------------------------------------------------
// A new material point class is defined to store the elastic parameters for 
// each integration point.
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt);

	//! The init function is used to intialize data
	void Init();

	//! copy material point data (for running restarts) \todo Is this still used?
	FEMaterialPoint* Copy();

	//! copy material point data (for running restarts) \todo Is this still used?
	void Serialize(DumpStream& dmp);

public:
	// These are the material parameters
	double		m_D;		// collagen density (?)
	double		m_DA;		// degree of anisotropy (?)

public:
	GridPoint	m_pt;	// grid point location of this material point

public:
	double vessel_weight;
	double matrix_weight;
	FEMaterialPoint* vessPt;
	FEMaterialPoint* matPt;

	DECLARE_PARAMETER_LIST();

public:
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};

//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterial : public FEElasticMaterial
{
public:
	struct SPROUT
	{
		vec3d		sprout;	// sprout direction
		FEElement*	pel;	// element in which this sprout lies
		double		r[3];	// iso-parameteric elements
	};

public:
	FEAngioMaterial(FEModel* pfem);

public:
	// material initialization
	bool Init();

	// Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp);

	// Calculate the active Angio stress
	mat3ds AngioStress(FEAngioMaterialPoint& mp);

	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData();

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp);

	double StrainEnergyDensity(FEMaterialPoint& mp);

public:
	// clear all sprouts
	void ClearSprouts();

	// add a sprout force
	// at position r with directional vector n
	void AddSprout(const vec3d& r, const vec3d& n);

	// return number of sprouts
	int Sprouts() { return (int) m_spr.size(); }

	// get a sprout
	SPROUT& GetSprout(int i) { return m_spr[i]; }

	// calculate the current spatial position, given an element and local coordinates
	vec3d CurrentPosition(FEElement* pe, double r, double s, double t);

	// we use this to define a sprout in the material section of the input file
	virtual void SetParameter(FEParam& p);

	//! Assign a grid
	void SetFEAngio(FEAngio* pangio) { m_pangio = pangio; }

private:
	double	m_a;
	double	m_b;
	double  m_N;
	int mat_id;

	// user-defined sprouts
	vec3d	m_s;	//!< dummy parameter used for reading sprouts from the input file
	vector<vec3d>	m_suser;

	vector<SPROUT>	m_spr;

	DECLARE_PARAMETER_LIST();

public:
	double scale;

public:
	int sym_planes[7];

	double Sx;
	double Sy;
	double Sz; 
	
	bool sym_on;

	vec3d sym;
	double sym_vects[7][3];

private:
	FEAngio * m_pangio;

public:
	void ApplySym();
	void MirrorSym(vec3d x, mat3ds &si, SPROUT sp, double den_scale);

private:
	FEPropertyT<FESolidMaterial> vessel_material;
	FEPropertyT<FESolidMaterial> matrix_material;
};

//-----------------------------------------------------------------------------
class FEPressureMaterial : public FEElasticMaterial
{
public:
	FEPressureMaterial(FEModel* pfem) : FEElasticMaterial(pfem){}

	mat3ds Stress(FEMaterialPoint& mp);

	tens4ds Tangent(FEMaterialPoint& mp);

public:
	double	m_p;	// pressure

	DECLARE_PARAMETER_LIST();
};