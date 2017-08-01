///////////////////////////////////////////////////////////////////////
// FESproutBodyForce.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FESPROUTBODYFORCE class is used to represent cell-generated forces 
// applied to the ECM
///////////////////////////////////////////////////////////////////////
#pragma once
#include "FEBioMech/FEBodyForce.h"
#include "FECore/FEParameterList.h"
#include "FECore/FEElement.h"


//REFACTOR: make some fields private
//-----------------------------------------------------------------------------
class FESproutBodyForce : public FEBodyForce
{
// Public field: SPROUT class
public:
	struct SPROUT
	{
		SPROUT(vec3d & r,vec3d & sp):rc(r), sprout(sp){}
		vec3d	rc;			// center of sprout force
		vec3d	sprout;		// sprout direction
	};

// Public functions:
	explicit FESproutBodyForce(FEModel* pfem);

	vec3d force(FEMaterialPoint& mp) override;
	mat3ds stiffness(FEMaterialPoint& mp) override;

	void Serialize(DumpStream& ar) override;

	void ApplySym();
	void MirrorSym(vec3d x, vec3d &f, SPROUT sp);

	bool Init() override;
	void Update() override;

	void AddSprout(vec3d r, vec3d s)
	{
		m_sp.emplace_back(r,s);
	}

	void ClearSprouts() { m_sp.clear(); }

	int Sprouts() const { return static_cast<int>(m_sp.size()); }

	SPROUT& GetSprout(int i) { return m_sp[i]; }

// Public fields:
	double	m_a, m_b;

	vector<SPROUT>	m_sp;

	double     m_factor;
	
	int		m_inode;

	bool	m_brigid;

	FESolidElement* m_pel;		//!< element in which point m_r0 lies
	double m_rs[3];				//!< isoparametric coordinates
		
	int sym_planes[7];

	double Sx;
	double Sy;
	double Sz; 
	
	bool sym_on;

	DECLARE_PARAMETER_LIST();

	vec3d sym;
	double sym_vects[7][3];
};
