#pragma once
#include <FECore/vec3d.h>
#include "FECore/FEElement.h"
#include "FECore/FESolidDomain.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "Segment.h"
#include "CultureParameters.h"
#include "ECMInitializer.h"
#include "Culture.h"
#include <FECore/FENormalProjection.h>
#include "FiberManager.h"


class FEAngioMaterialPoint;

//contains the shared functionality, material pointers may need to be passed in to get portions of the functionality
class FEAngioMaterialBase
{
public:
	struct SPROUT
	{
		explicit SPROUT(const vec3d & dir, FESolidElement * el, double * local, FEAngioMaterialBase * m0, FEElasticMaterial *m1);
		vec3d		sprout;	// sprout direction
		FESolidElement*	pel;	// element in which this sprout lies
		double		r[3];	// iso-parameteric elements
		FEAngioMaterialBase * mat0;
		FEElasticMaterial * mat1;
	};
	FEAngioMaterialBase();
	virtual ~FEAngioMaterialBase(){}

	//bool SharedInit();

	// return number of sprouts
	int Sprouts() const { return (int)m_spr.size(); }

	// get a sprout
	SPROUT& GetSprout(int i) { return m_spr[i]; }

	void CreateSprouts(double scale, FEElasticMaterial* emat);

	void UpdateSprouts(double scale, FEElasticMaterial* emat);

	// clear all sprouts
	void ClearSprouts();

	// add a sprout force
	// at position r with directional vector n
	void AddSprout(const vec3d& r, const vec3d& n, FESolidDomain * domain, int elemindex, FEElasticMaterial* emat);
	void AddSprout(const vec3d& r, const vec3d& n, FEDomain * domain, FEElasticMaterial* emat);
	void AddSprout(const Segment::TIP & tip, FEElasticMaterial* emat);


	bool FindGridPoint(const vec3d & r, GridPoint & p) const;

	bool FindGridPoint(const vec3d & r, FESolidDomain * domain, int elemindex, GridPoint & p) const;

	// calculate the current spatial position, given an element and local coordinates
	vec3d CurrentPosition(FESolidElement * pe, double r, double s, double t) const;

	void AdjustMeshStiffness(FEMaterial* mat);

	void UpdateFiberManager();

	//! Assign a grid
	void SetFEAngio(FEAngio* pangio) { m_pangio = pangio; }

	void MirrorSym(vec3d x, mat3ds &si, SPROUT sp, double den_scale);

	void UpdateSproutStressScaling();

	bool InitCulture();

	void Update();

	bool Overwrite() const;

	double GetAnisotropy() const;

	//begin virtual functions
	virtual void InitializeFibers()=0;

	virtual mat3ds AngioStress(FEAngioMaterialPoint& mp)=0;

	virtual void FinalizeInit()=0;

	virtual void UpdateECM()=0;

	virtual void UpdateGDMs()=0;

	virtual bool InitECMDensity(FEAngio * angio)=0;

	FEAngio * m_pangio;
	CultureParameters m_cultureParams;
	Culture * m_cult;

	ECMInitializer * ecm_initializer;

	double scale;

	int sym_planes[7];

	bool sym_on;

	vec3d sym;
	double sym_vects[7][3];

	std::vector<int> domains;
	std::vector<FEDomain*> domainptrs;
	std::unordered_map<int, int> node_map;//maps global node number to domain local node numbers
	std::unordered_map<FEDomain *, int> meshOffsets;
	FESurface * exterior_surface;
	FENormalProjection * normal_proj;

	FiberManager * fiber_manager;

	int mat_id;

	// user-defined sprouts
	vec3d	m_s;	//!< dummy parameter used for reading sprouts from the input file
	std::vector<vec3d>	m_suser;

	//m_spr is the underlying storage for sprouts
	std::vector<SPROUT>	m_spr;
	KDTree<std::pair<size_t, std::vector<SPROUT> *>, std::vector<double>> sprouts;
protected:
};