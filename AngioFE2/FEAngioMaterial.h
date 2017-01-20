#pragma once
#include "StdAfx.h"
#include <FEBioMech/FEElasticMaterial.h>
#include <FECore/FEDataArray.h>
#include <FECore/FESurface.h>
#include <FECore/FENormalProjection.h>
#include "FEAngio.h"
#include "Culture.h"
#include "FEProbabilityDistribution.h"

//-----------------------------------------------------------------------------
// A new material point class is defined to store the elastic parameters for 
// each integration point.
class FEAngioMaterialPoint : public FEMaterialPoint
{
public:
	//! constructor
	FEAngioMaterialPoint(FEMaterialPoint* pt, FEMaterialPoint* vesselPt, FEMaterialPoint *matrixPt);

	//! The init function is used to intialize data
	void Init() override;

	//! copy material point data (for running restarts) todo Is this still used?
	FEMaterialPoint* Copy() override;

	//! copy material point data (for running restarts) todo Is this still used?
	void Serialize(DumpStream& dmp) override;

	// These are the material parameters
	double		m_D;		// collagen density (?)
	double		m_DA;		// degree of anisotropy (?)
	//TODO: consider adding a custom weigth parameter per element 

	GridPoint	m_pt;	// grid point location of this material point

	double vessel_weight;
	double matrix_weight;
	FEMaterialPoint* vessPt;
	FEMaterialPoint* matPt;

	DECLARE_PARAMETER_LIST();

public:
	static FEAngioMaterialPoint* FindAngioMaterialPoint(FEMaterialPoint* mp);
};
class FEAngioMaterial;
//the interface for the seeding of the extra cellular matrix
//interfaces can be used for a number of objectives
//these include controlling what happens on the boundaries and seeding according to a distribution

//needs to set the intitial ecm density and anisotropy
//if the density has not yet been modified it will be 0.0
class ECMInitializer
{
public:
	virtual ~ECMInitializer(){}
	virtual void seedECMDensity(FEAngioMaterial * mat) = 0;
	virtual bool overwrite(){ return true; }
	virtual void updateECMdensity(FEAngioMaterial * mat);
};
class ECMInitializerConstant : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
};
class ECMInitializerSpecified : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
};
class ECMInitializerNoOverwrite : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
	bool overwrite() override{ return false; }
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
	FEAngioMaterial(FEModel* pfem);
	virtual ~FEAngioMaterial();

	friend class Fileout;
	friend class FEPlotMatrixStress;
	friend class FEPlotVesselStress;
	friend class FEPlotMatrixTangent;

	// material initialization
	bool Init() override;
	void FinalizeInit();

	void AdjustMeshStiffness();

	bool InitCollagenFibers();

	void CreateSprouts(double scale);

	void UpdateSprouts(double scale);

	void UpdateSproutStressScaling();

	bool InitCulture();

	void Update();

	void UpdateECM();

	bool Overwrite() const;

	vec3d CollagenDirection(GridPoint& pt);

	bool InitECMDensity(FEAngio * angio);

	// Calculate Cauchy-stress
	mat3ds Stress(FEMaterialPoint& mp) override;

	// Calculate the active Angio stress
	mat3ds AngioStress(FEAngioMaterialPoint& mp);

	// Calculate spatial elasticity tangent
	tens4ds Tangent(FEMaterialPoint& mp) override;

	//! create material point data for this material
	FEMaterialPoint* CreateMaterialPointData() override;

	//! Set the local coordinate system for a material point (overridden from FEMaterial)
	void SetLocalCoordinateSystem(FEElement& el, int n, FEMaterialPoint& mp) override;

	double StrainEnergyDensity(FEMaterialPoint& mp) override;

	// clear all sprouts
	void ClearSprouts();

	// add a sprout force
	// at position r with directional vector n
	void AddSprout(const vec3d& r, const vec3d& n, FEDomain * domain, int elemindex);
	void AddSprout(const vec3d& r, const vec3d& n, FEDomain * domain);
	void AddSprout(const Segment::TIP & tip);

	// return number of sprouts
	int Sprouts() const { return (int) m_spr.size(); }

	// get a sprout
	SPROUT& GetSprout(int i) { return m_spr[i]; }

	// calculate the current spatial position, given an element and local coordinates
	vec3d CurrentPosition(FEElement* pe, double r, double s, double t) const;

	// we use this to define a sprout in the material section of the input file
	void SetParameter(FEParam& p) override;

	//! Assign a grid
	void SetFEAngio(FEAngio* pangio) { m_pangio = pangio; }

	double GetAnisotropy() const;
	
	void SetupSurface();

	bool FindGridPoint(const vec3d & r, GridPoint & p) const;
	bool FindGridPoint(const vec3d & r, FEDomain * domain, int elemindex, GridPoint & p) const;
private:
	ECMInitializer * ecm_initializer;

	

	int mat_id;

	// user-defined sprouts
	vec3d	m_s;	//!< dummy parameter used for reading sprouts from the input file
	vector<vec3d>	m_suser;

	vector<SPROUT>	m_spr;

	DECLARE_PARAMETER_LIST();

public:
	double scale;

	int sym_planes[7];
	
	bool sym_on;

	vec3d sym;
	double sym_vects[7][3];

	std::vector<int> domains;
	std::vector<FEDomain*> domainptrs;
	std::unordered_map<FEDomain *,int> meshOffsets;
	FESurface * exterior_surface;
	FENormalProjection * normal_proj;
	Culture * m_cult;
	FEAngio * m_pangio;
	CultureParameters m_cultureParams;
private:
	

	

public:
	void ApplySym();
	void MirrorSym(vec3d x, mat3ds &si, SPROUT sp, double den_scale);

private:
	FEPropertyT<FESolidMaterial> vessel_material;
	FEPropertyT<FESolidMaterial> matrix_material;
	FEPropertyT<FEProbabilityDistribution> length_to_branch;
	FEPropertyT<FragmentBranching> fbrancher;
};

//-----------------------------------------------------------------------------
class FEPressureMaterial : public FEElasticMaterial
{
public:
	FEPressureMaterial(FEModel* pfem) : FEElasticMaterial(pfem){ m_p = 0; }

	mat3ds Stress(FEMaterialPoint& mp) override;

	tens4ds Tangent(FEMaterialPoint& mp) override;

	double	m_p;	// pressure

	DECLARE_PARAMETER_LIST();
};