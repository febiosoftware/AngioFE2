#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <cmath>
#include "Elem.h"
#include "Culture.h"

//-----------------------------------------------------------------------------
// define the material parameters
BEGIN_PARAMETER_LIST(FEAngioMaterialPoint, FEMaterialPoint)
	ADD_PARAMETER(m_D, FE_PARAM_DOUBLE, "dens");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterialPoint::FEAngioMaterialPoint(FEMaterialPoint* pt) : FEMaterialPoint(pt)
{
	m_D = 0.0;
}

//-----------------------------------------------------------------------------
//! The init function is used to intialize data
void FEAngioMaterialPoint::Init(bool bflag)
{
}

//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
FEMaterialPoint* FEAngioMaterialPoint::Copy()
{
	FEAngioMaterialPoint* pt = new FEAngioMaterialPoint(*this);
	if (m_pNext) pt->m_pNext = m_pNext->Copy();
	return pt;
}

//-----------------------------------------------------------------------------
//! copy material point data (for running restarts) \todo Is this still used?
void FEAngioMaterialPoint::Serialize(DumpStream& dmp)
{
	if (dmp.IsSaving())
	{
		dmp << m_D;
	}
	else
	{
		dmp >> m_D;
	}
	FEMaterialPoint::Serialize(dmp);
}

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_a, FE_PARAM_DOUBLE, "a");
	ADD_PARAMETER(m_b, FE_PARAM_DOUBLE, "b");
	ADD_PARAMETER(m_N, FE_PARAM_DOUBLE, "N");
	ADD_PARAMETER(m_s, FE_PARAM_VEC3D , "sprout");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	scale = 1.0;

	m_pangio = 0;

	// Create symmetry vectors
	sym_planes[0] = 0; sym_planes[1] = 0; sym_planes[2] = 0; sym_planes[3] = 0; sym_planes[4] = 0; sym_planes[5] = 0; sym_planes[6] = 0;
	Sx = 0.; Sy = 0.; Sz = 0.;
	
	sym_vects[0][0] = 1.; sym_vects[0][1] = 0.; sym_vects[0][2] = 0.;
	sym_vects[1][0] = 0.; sym_vects[1][1] = 1.; sym_vects[1][2] = 0.;
	sym_vects[2][0] = 0.; sym_vects[2][1] = 0.; sym_vects[2][2] = 1.;
	sym_vects[3][0] = 1.; sym_vects[3][1] = 1.; sym_vects[3][2] = 0.;
	sym_vects[4][0] = 1.; sym_vects[4][1] = 0.; sym_vects[4][2] = 1.;
	sym_vects[5][0] = 0.; sym_vects[5][1] = 1.; sym_vects[5][2] = 1.;
	sym_vects[6][0] = 1.; sym_vects[6][1] = 1.; sym_vects[6][2] = 1.;

	sym_on = false;
}

//-----------------------------------------------------------------------------
bool FEAngioMaterial::Init()
{
	if (FEElasticMaterial::Init() == false) return false;

	// add the user sprouts
	for (int i=0; i<m_suser.size(); ++i)
	{
		AddSprout(m_suser[i], vec3d(0,0,0));
	}
	m_suser.clear();

	// initialize material point data
	vec3d x[FEElement::MAX_NODES];
	FEMesh& mesh = GetFEModel()->GetMesh();
	Grid& grid = m_pangio->GetGrid();
	for (int n=0; n<mesh.Domains(); ++n)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(n));
		FEMaterial* pm = dom.GetMaterial();
		FEAngioMaterial* pam = dynamic_cast<FEAngioMaterial*>(pm->FindComponentByType("angio"));
		if (pam == this)
		{
			// loop over all elements
			int NE = dom.Elements();
			for (int i=0; i<NE; ++i)
			{
				// get the next element
				FEElement& el = dom.Element(i);
				int neln = el.Nodes();
				
				// get the nodal coordinates
				for (int j=0; j<neln; ++j) x[j] = mesh.Node(el.m_node[j]).m_rt;

				// loop over all integration points
				int nint = el.GaussPoints();
				for (int j=0; j<nint; ++j)
				{
					FEMaterialPoint& mp = *el.GetMaterialPoint(j);

					int nc = mp.Components();
					for (int k=0; k<nc; ++k)
					{
						FEMaterialPoint& mk = *mp.GetPointData(k);
					
						FEAngioMaterialPoint* ap = mk.ExtractData<FEAngioMaterialPoint>();
						if (ap)
						{
							// calculate the spatial location of the integration point
							vec3d r = el.Evaluate(x, j);

							// calculate the GridPoint data for this point.
							if (grid.FindGridPoint(r, ap->m_pt) == false)
							{
								return false;
							}
						}
					}
				}
			}
		}
	}

	return true;
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.m_szname, "sprout") == 0)
	{
		m_suser.push_back(m_s);
	}
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::ClearSprouts()
{
	m_spr.clear();
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));

	SPROUT s;
	s.pel = dom.FindElement(r, s.r);
	s.sprout = t;

	assert(s.pel);

	m_spr.push_back(s);
}

//-----------------------------------------------------------------------------
vec3d FEAngioMaterial::CurrentPosition(FEElement* pe, double r, double s, double t)
{
	double H[8];
	H[0] = 0.125*(1.0 - r)*(1.0 - s)*(1.0 - t);
	H[1] = 0.125*(1.0 + r)*(1.0 - s)*(1.0 - t);
	H[2] = 0.125*(1.0 + r)*(1.0 + s)*(1.0 - t);
	H[3] = 0.125*(1.0 - r)*(1.0 + s)*(1.0 - t);
	H[4] = 0.125*(1.0 - r)*(1.0 - s)*(1.0 + t);
	H[5] = 0.125*(1.0 + r)*(1.0 - s)*(1.0 + t);
	H[6] = 0.125*(1.0 + r)*(1.0 + s)*(1.0 + t);
	H[7] = 0.125*(1.0 - r)*(1.0 + s)*(1.0 + t);

	FEMesh& mesh = GetFEModel()->GetMesh();
	vec3d rc(0,0,0);
	for (int i=0; i<8; ++i) rc += mesh.Node(pe->m_node[i]).m_rt*H[i];

	return rc;
}

//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& ep = *mp.ExtractData<FEElasticMaterialPoint>();
	FEAngioMaterialPoint& ap = *mp.ExtractData<FEAngioMaterialPoint>();

	// get density scale factor
	Culture& cult = m_pangio->GetCulture();
	double den_scale = cult.FindDensityScale(ap.m_pt);

	// loop over all sprout tips
	mat3ds s(0.0);
	int NS = Sprouts();

//#pragma omp parallel for shared(s)
	for (int i=0; i<NS; ++i)
	{
		SPROUT& sp = m_spr[i];
		
		// current position of sprout force
		vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);

		// current position of integration point
		vec3d y = ep.m_rt;

		vec3d r = y - x;
		double l = r.unit();

		sp.sprout.unit();															// Normalize the sprout direction vector
			
		double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

		double p = den_scale*scale*m_a*(pow(cos(theta/2),m_N))*exp(-m_b*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

		//double p = m_a*exp(-m_b*l);

		mat3ds si = dyad(r)*p;

		if (sym_on == true)															// If symmetry is turned on, apply symmetry
			MirrorSym(y, si, sp, den_scale);

//#pragma omp critical
		s += si;
	}
	return s;
}

//-----------------------------------------------------------------------------
tens4ds FEAngioMaterial::Tangent(FEMaterialPoint& mp)
{
	tens4ds C(0.0);
	return C;
}

//=============================================================================
BEGIN_PARAMETER_LIST(FEPressureMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_p, FE_PARAM_DOUBLE, "p");
END_PARAMETER_LIST();

mat3ds FEPressureMaterial::Stress(FEMaterialPoint& pt)
{
	mat3dd I(1.0);
	return I*m_p;
}

tens4ds FEPressureMaterial::Tangent(FEMaterialPoint& pt)
{
	return tens4ds(0.0);
}

//============================================================================

///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - ApplySym
//      Determine if symmetry is turned on, if so create the symmetry vectors
///////////////////////////////////////////////////////////////////////

void FEAngioMaterial::ApplySym()
{
	if (Sx != 0)															// Turn on x symmetry
		sym_planes[0] = 1;

	if (Sy != 0)															// Turn on y symmetry
		sym_planes[1] = 1;

	if (Sz != 0)															// Turn on z symmetry
		sym_planes[2] = 1;

	if ((Sx != 0) && (Sy != 0))												// Turn on x and y symmetry
		sym_planes[3] = 1;

	if ((Sx != 0) && (Sz != 0))												// Turn on x and z symmetry
		sym_planes[4] = 1;

	if ((Sy != 0) && (Sz != 0))												// Turn on y and z symmetry
		sym_planes[5] = 1;

	if ((Sx != 0) && (Sy != 0) && (Sz != 0))								// Turn on x y and z symmetry
		sym_planes[6] = 1;
	
	if (sym_planes[0] + sym_planes[1] + sym_planes[2] + sym_planes[3] + sym_planes[4] + sym_planes[5] + sym_planes[6] != 0)
		sym_on = true;														

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - MirrorSym
//      Calculate force due to mirrored vessels at a particular material point at position x
///////////////////////////////////////////////////////////////////////

void FEAngioMaterial::MirrorSym(vec3d y, mat3ds &s, SPROUT sp, double den_scale)
{
	sym.x = Sx; sym.y = Sy; sym.z = Sz;										// Set the position of the symmetry planes
	mat3ds ssym; ssym.zero();
	vec3d sprout_vect;														// Sprout direction vector
	vec3d sym_v;															// Symmetry vector

	for (int i = 0; i < 7; i++){											// For each of the possible symmetry planes
		if (sym_planes[i] == 1){												// If that symmetry plane is turned on...
			sym_v.x = sym_vects[i][0]; sym_v.y = sym_vects[i][1]; sym_v.z = sym_vects[i][2];	// Obtain the symmetry vector

			vec3d x = CurrentPosition(sp.pel, sp.r[0], sp.r[1], sp.r[2]);
			vec3d r = y - x;																	// Draw the vector r from the sprout location to the position of the material point
			
			//r.x += 2*sym_v.x*(sym.x - x.x); r.y += 2*sym_v.y*(sym.y - x.y); r.z += 2*sym_v.z*(sym.z - x.z);		// Find r for the mirrored vessel sprout
			r.x = r.x + sym_v.x*sym.x; r.y = r.y + sym_v.y*sym.y; r.z = r.z + sym_v.z*sym.z;
			double l = r.unit();													// Find the length of r
			
			sprout_vect.x = sp.sprout.x; sprout_vect.y = sp.sprout.y; sprout_vect.z = sp.sprout.z;	// Set the sprout direction vector 
			sprout_vect.unit();														// Normalize the sprout direction vector
			
			if (m_N != 0){														// If a directional sprout force is being used...
				switch (i) {
				case 0:
					sprout_vect.x = -sprout_vect.x; break;								// Mirror across x
				case 1:
					sprout_vect.y = -sprout_vect.y; break;								// Mirror across y
				case 2:
					sprout_vect.z = -sprout_vect.z; break;								// Mirror across z
				case 3:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; break;	// Mirror across x and y
				case 4:
					sprout_vect.x = -sprout_vect.x; sprout_vect.z = -sprout_vect.z; break;	// Mirror across x and z
				case 5:
					sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;	// Mirror across y and z
				case 6:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;}	// Mirror across x y and z
			}
	
			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double p = den_scale*scale*m_a*(pow(cos(theta/2),m_N))*exp(-m_b*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation
			
			if ((p != p) || (r.x != r.x) || (r.y != r.y) || (r.z != r.z)){			// If the mirrored force vector isn't real...
				p = 0.; r.x = 0.; r.y = 0.; r.z = 0.;}									// Set it to zero

			ssym += dyad(r)*p;
		}
	}


	s += ssym;															// Add the symmetry results to the force vector

	return;
}

//-----------------------------------------------------------------------------
FEMaterialPoint* FEAngioMaterial::CreateMaterialPointData()
{
	return new FEAngioMaterialPoint(new FEElasticMaterialPoint); 
}
