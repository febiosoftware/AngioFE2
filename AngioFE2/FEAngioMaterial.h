#pragma once
#include <FEBioMech/FEElasticMaterial.h>
#include "Grid.h"

//-----------------------------------------------------------------------------
// A new material point class is defined to store the elastic parameters for 
// each integration point.
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEAngioMaterialPoint(FEMaterialPoint* pt);

	//! The init function is used to intialize data
	virtual void Init(bool bflag);

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual FEMaterialPoint* Copy();

	//! copy material point data (for running restarts) \todo Is this still used?
	virtual void ShallowCopy(DumpStream& dmp, bool bsave);

public:
	// These are the material parameters
	double	m_D;

	DECLARE_PARAMETER_LIST();
};

//-----------------------------------------------------------------------------
// Class implementing a stress induced by a non-local point force
class FEAngioMaterial : public FEElasticMaterial
{
public:
	struct SPROUT
	{
		bool	bactive;
		vec3d	sprout;

		FEElement*	pel;
		double	r[3];	// iso-parameteric elements
	};

public:
	FEAngioMaterial(FEModel* pfem);

public:
	// material initialization
	bool Init();

	// Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp);

	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp);

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData();

public:
	// add a sprout force
	void AddSprout(const vec3d& r, const vec3d& n);

	// update the sprout position
	void UpdateSprout(SPROUT& s, const vec3d& r);

	// get the total sprout stress
	mat3ds SproutStress(const vec3d& x);

	// return number of sprouts
	int Sprouts() { return (int) m_spr.size(); }

	// get a sprout
	SPROUT& GetSprout(int i) { return m_spr[i]; }

	vec3d CurrentSproutPosition(SPROUT& s);

	// we use this to define a sprout in the material section of the input file
	virtual void SetParameter(FEParam& p);

	//! Assign a grid (TODO: Do we really need this?)
	void SetGrid(Grid* pgrid) { m_pgrid = pgrid; }

private:
	double	m_a;
	double	m_b;
	double  m_N;

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
	Grid * m_pgrid;

public:
	void ApplySym();
	void MirrorSym(vec3d x, mat3ds &si, SPROUT sp, double den_scale);
	double findDenScale(double xpt, double ypt, double zpt);
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
