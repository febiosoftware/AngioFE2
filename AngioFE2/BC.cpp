#include "StdAfx.h"
#include "BC.h"
#include "Segment.h"
#include "FEAngio.h"
#include "FECore/FEMesh.h"
#include "FECore/FEModel.h"
#include "Culture.h"
#include <cassert>
#include <cmath>
#include "angio3d.h"

//-----------------------------------------------------------------------------
BC::BC(FEModel * model) : FEMaterial(model)
{
	AddProperty(&mbc, "mbc");
}  

void BC::SetCulture(Culture * cp)
{
	culture = cp;
	mbc->SetCulture(cp);
}
void MBC::SetCulture(Culture * cp)
{
	culture = cp;
}

//-----------------------------------------------------------------------------
BC::~BC()
{
}

//-----------------------------------------------------------------------------
// Checks if a newly-created segment violates the boundary faces of the element in which it occupies
void BC::CheckBC(Segment &seg)
{
	
	{
		auto rseg = culture->RecentSegments();
		if(rseg.size() >= culture->m_cultParams->max_recursion_depth)
		{
			return;
		}
	}
	//new implementation may be run on all segments will add the segment once the boundaries are safe
	//remove the check of whether or not to check the boundary condition in Culture
	assert(seg.tip_c(0).pt.nelem != -1);
	//make sure that the defining tip has not drifted too far
	vec3d pq = culture->m_pmat->m_pangio->Position(seg.tip_c(0).pt);
	assert((culture->m_pmat->m_pangio->Position(seg.tip_c(0).pt) - seg.tip_c(0).pt.r).norm() < 1.0);
	//if the second tip is unititialized make sure it is not in the same element as tip(0)
	FESolidElement * tse = &seg.tip_c(0).pt.ndomain->Element(seg.tip_c(0).pt.elemindex);

	if (seg.tip_c(1).pt.nelem == -1)
	{
		
		double arr[3];
		if (tse && culture->m_pmat->m_pangio->IsInsideHex8(tse, seg.tip_c(1).pt.r, culture->m_pmat->m_pangio->GetMesh(), arr))
		{
			//just copy the data from tip(0)
			seg.tip(1).pt.q.x = arr[0]; seg.tip(1).pt.q.y = arr[1]; seg.tip(1).pt.q.z = arr[2];
			seg.tip(1).pt.nelem = seg.tip_c(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip_c(0).pt.ndomain;
			seg.tip(1).pt.elemindex = seg.tip_c(0).pt.elemindex;
		}
	}
	//this is the costliest part of the boundary check if reached
	if (seg.tip_c(1).pt.nelem == -1)
	{
		culture->m_pmat->FindGridPoint(seg.tip_c(1).pt.r ,seg.tip(1).pt);
	}
	//check if the segment is in another angio material
	if (seg.tip_c(1).pt.nelem == -1)
	{
		FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
		double r[3];
		FESolidElement * se = mesh->FindSolidElement(seg.tip_c(1).pt.r, r);
		
		if (se)
		{
			FEAngioMaterialBase * angm = culture->m_pmat->m_pangio->GetAngioComponent(culture->m_pmat->m_pangio->m_fem->GetMaterial(se->GetMatID()));
			if(angm)
			{
				GridPoint & cpt = seg.tip(1).pt;
				cpt.q = vec3d(r[0], r[1], r[2]);
				cpt.ndomain = dynamic_cast<FESolidDomain*>(se->GetDomain());
				cpt.nelem = se->GetID();
				cpt.elemindex = se->GetID() - 1 - angm->meshOffsets.at(cpt.ndomain);
				seg.Update();
				if (mbc->acceptBoundary(culture->m_pmat, angm) && (angm != this->culture->m_pmat))
				{
					angm->m_cult->ClearRecents();
					mbc->handleBoundary(culture->m_pmat, angm, seg);
					auto orseg = angm->m_cult->RecentSegments();
					for (size_t i = 0; i < orseg.size(); i++)
					{
						culture->AddToRecents(orseg[i]);
					}
					return;
				}
			}
		}
	}
	
	//if both are in the same element just add the segment
	if (seg.tip_c(0).pt.nelem == seg.tip_c(1).pt.nelem)
	{
		seg.Update();
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			culture->AddSegment(seg);
			return;
		}
	}
	//if both elmeents are in the same material add the segment
	//this implies that angio materials are convex or that the user will in some way mitigate the nonconvex portions of the material
	if (seg.tip_c(1).pt.ndomain != nullptr && seg.tip_c(1).pt.elemindex >= 0)
	{
		FEElement & t1se = seg.tip_c(1).pt.ndomain->ElementRef(seg.tip_c(1).pt.elemindex);
		FEElement & t0se = seg.tip_c(0).pt.ndomain->ElementRef(seg.tip_c(0).pt.elemindex);
		if (t1se.GetMatID() == t0se.GetMatID())
		{
			if (seg.length() >= culture->m_cultParams->min_segment_length)
			{
				culture->AddSegment(seg);
				return;
			}
		}
	}
	
	vec3d dir = seg.tip_c(1).pt.r - seg.tip_c(0).pt.r;
	double dist = dir.norm();
	assert(dir.norm() > 0.0);

	dir.unit();

	FESolidDomain & sd = reinterpret_cast<FESolidDomain&>(*seg.tip_c(0).pt.ndomain);
	FESolidElement * se = &sd.Element(seg.tip_c(0).pt.elemindex);
	assert(se);//make sure we have an element

	FEMaterial * mat = sd.GetMaterial();
	vec3d lastgood_pt;
	double rs[2];
	double g;
	assert(culture->m_pmat);
	FESurface * surf = culture->m_pmat->exterior_surface;
	if (culture->m_pmat->m_pangio->m_fe_element_data[se->GetID()].surfacesIndices.size() == 0)
		return;
	int hcount = 0;//counts how many times the boundary is handled
	std::vector<int> & edinices = culture->m_pmat->m_pangio->m_fe_element_data[se->GetID()].surfacesIndices;
	culture->m_pmat->normal_proj->SetSearchRadius(seg.length() * 2 + culture->m_cultParams->min_segment_length*2);
	for (size_t i = 0; i < edinices.size(); i++)
	{
		if (surf->Elements() > culture->m_pmat->m_pangio->m_fe_element_data[se->GetID()].surfacesIndices[i])
		{
			FESurfaceElement & surfe = reinterpret_cast<FESurfaceElement&>(surf->ElementRef(culture->m_pmat->m_pangio->m_fe_element_data[se->GetID()].surfacesIndices[i]));
			if (culture->m_pmat->exterior_surface->Intersect(surfe, seg.tip_c(0).pt.r, dir, rs, g, 0.0001))//see the epsilon in FIndSolidElement
			{
				lastgood_pt = surf->Local2Global(surfe, rs[0], rs[1]);
				double rdist = (lastgood_pt - seg.tip_c(1).pt.r).norm();
				if((abs(g) <= seg.expected_length) && (rdist < 1.01*(seg.expected_length)))
				{
					//set last_goodpt
					
					vec3d projected_pt = (seg.tip_c(1).pt.r - seg.tip_c(0).pt.r);
					projected_pt.unit();
					projected_pt = seg.tip_c(0).pt.r + projected_pt*g;

					seg.tip(1).pt.elemindex = seg.tip_c(0).pt.elemindex;
					seg.tip(1).pt.nelem = seg.tip_c(0).pt.nelem;
					seg.tip(1).pt.ndomain = seg.tip_c(0).pt.ndomain;
					return HandleBoundary(seg, lastgood_pt, rs, se);
				}
				
			}
		}
		
		
	}
	FESurfaceElement * surfe = culture->m_pmat->normal_proj->Project(seg.tip_c(0).pt.r, dir, rs);
	if (!surfe)
	{
		//printf("no surface element found\n");
		//assert(false);
		/*
		surfe = culture->m_pmat->normal_proj->Project3(seg.tip(0).pt.r, -dir, rs);
		if (surfe)
		{
			lastgood_pt = surf->Local2Global(*surfe, rs[0], rs[1]);

			seg.tip(1).pt.elemindex = seg.tip(0).pt.elemindex;
			seg.tip(1).pt.nelem = seg.tip(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			return HandleBoundary(seg, lastgood_pt, rs, se);
		}
		else
		{
			printf("segment adjustment failed\n");
		}
		*/
	}
	else
	{

		int eindex = surfe->m_elem[0];
			
		vec3d pos = surf->Local2Global(*surfe, rs[0], rs[1]);
		vec3d dir = pos - seg.tip(0).pos();
		dir.unit();

		seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
		seg.tip(1).pt.elemindex = eindex - culture->m_pmat->meshOffsets.at(seg.tip(1).pt.ndomain);
		seg.tip(1).pt.nelem = eindex + 1;
		seg.tip(1).pt.r = seg.tip(0).pt.r + dir*seg.length();
		FESolidElement & se = reinterpret_cast<FESolidElement&>(seg.tip(1).pt.ndomain->ElementRef(seg.tip(1).pt.elemindex));
		
		return HandleBoundary(seg,  pos, rs, &se);
	}
} 


bool BC::ChangeOfMaterial(Segment & seg) const
{
	//if either is not in the model they are in different materials
	if (seg.tip(0).pt.nelem == -1 || seg.tip(1).pt.nelem == -1)
		return true;
	//if both are in the same element return
	if (seg.tip(0).pt.nelem == seg.tip(1).pt.nelem)
		return false;
	
	vec3d dir = seg.tip(1).pt.r - seg.tip(0).pt.r;

	//not sure a meaningful answer can be returned if the length is zero and is probably an error at the call site 
	assert(dir.norm() > 0.0);

	dir.unit();

	FESolidDomain & sd = reinterpret_cast<FESolidDomain&>(*seg.tip(0).pt.ndomain);
	FESolidElement * se = &sd.Element(seg.tip(0).pt.elemindex);
	assert(se);//make sure we have an element

	return false;
}

BEGIN_PARAMETER_LIST(BC, FEMaterial)
ADD_PARAMETER(angio_boundary_groups, FE_PARAM_INT, "angio_boundary_groups");
END_PARAMETER_LIST();
SameMBC::SameMBC(FEModel * model) : MBC(model)
{
	
}

StopBC::StopBC(FEModel * model) : BC(model)
{
	
}
void StopBC::HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se)
{
	//fill in the pt's data and add the segment
	//remaining distance is ignored
	//does not need to do anything for branching segments as the branch will 
	seg.tip(1).pt.r = lastGoodPt;
	seg.SetFlagOn(Segment::BC_DEAD);
	if (culture->m_pmat->FindGridPoint(lastGoodPt, seg.tip(1).pt.ndomain, seg.tip(1).pt.elemindex, seg.tip(1).pt))
	{
		assert(seg.tip_c(1).pt.nelem != -1);
		seg.Update();
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			seg.tip(0).bactive = false;
			seg.tip(1).bactive = false;
			if (seg.length() >= culture->m_cultParams->min_segment_length)
			{
				culture->AddSegment(seg);
				return;
			}
		}
	}

	else
	{
		//this denotes a mismatch between the FindSolidElement and FESurface output
		vec2d nrs(rs[0], rs[1]);
		seg.tip(1).pt.q = culture->m_pmat->m_pangio->FindRST(lastGoodPt, nrs, se);
		seg.tip(1).pt.r = lastGoodPt;
		seg.tip(1).pt.nelem = se->GetID();
		seg.tip(1).pt.ndomain = dynamic_cast<FESolidDomain*>(culture->m_pmat->domainptrs[0]); //hack
		seg.tip(1).pt.elemindex = se->GetID() - 1 - culture->m_pmat->meshOffsets.at(seg.tip(1).pt.ndomain);
		seg.tip(1).pt.r = culture->m_pmat->m_pangio->Position(seg.tip(1).pt);
		seg.Update();
		seg.tip(0).bactive = false;
		seg.tip(1).bactive = false;
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			culture->AddSegment(seg);
		}
	}
}
BouncyBC::BouncyBC(FEModel* model) : BC(model)
{
	
}


void BouncyBC::HandleBoundary(Segment & seg, vec3d lastGoodPt, double * rs, FESolidElement * se)
{
	//fill in the pt's data and add the segment
	//remaining distance is ignored
	double rdist = (lastGoodPt - seg.tip(1).pt.r).norm();
	assert(rdist < 1.01*(seg.expected_length));
	bool found0 = false;
	bool found1 = false;
	FEMesh * mesh = culture->m_pmat->m_pangio->GetMesh();
	assert(rdist >= 0.0);
	vec3d prevhead = seg.tip(1).pt.r;
	vec3d segadjdir = seg.tip(1).pt.r - seg.tip(0).pt.r;
	segadjdir.unit();
	seg.tip(1).pt.r = lastGoodPt;
	seg.SetFlagOn(Segment::BC_DEAD);
	if (culture->m_pmat->FindGridPoint(lastGoodPt, seg.tip(1).pt.ndomain, seg.tip(1).pt.elemindex, seg.tip(1).pt))
	{
		assert(seg.tip_c(1).pt.nelem != -1);
		seg.Update();
		if (seg.length() >= culture->m_cultParams->min_segment_length && seg.length() != 0.0)
		{
			seg.tip(0).bactive = false;
			seg.tip(1).bactive = false;
			culture->AddSegment(seg);
			found0 = true;
		}
	}

	else
	{
		//this denotes a mismatch between the FindSolidElement and FESurface output
		//printf("boundary segment out of bounds attempting to correct this\n");
		vec2d nrs(rs[0], rs[1]);
		seg.tip(1).pt.q = culture->m_pmat->m_pangio->FindRST(lastGoodPt, nrs, se);
		seg.tip(1).pt.r = lastGoodPt;
		seg.tip(1).pt.nelem = se->GetID();
		seg.tip(1).pt.ndomain = dynamic_cast<FESolidDomain*>(culture->m_pmat->domainptrs[0]); //hack
		seg.tip(1).pt.elemindex = se->GetID() - 1 - culture->m_pmat->meshOffsets.at(seg.tip(1).pt.ndomain);
		seg.tip(1).pt.r = culture->m_pmat->m_pangio->Position(seg.tip(1).pt);
		seg.Update();
		seg.tip(0).bactive = false;
		seg.tip(1).bactive = false;
		if (seg.length() >= culture->m_cultParams->min_segment_length && seg.length() != 0.0)
		{
			culture->AddSegment(seg);
			found0 = true;
		}
	}
	assert(seg.tip_c(1).pt.nelem != -1);
	vec3d dir = seg.tip(1).pt.r - seg.tip(0).pt.r;
	dir.unit();

	FESurfaceElement * surfe = culture->m_pmat->normal_proj->Project(seg.tip(0).pt.r, -dir, rs);
	if (surfe)
	{
		Segment reflSeg;
		reflSeg.tip(0).pt = seg.tip(1).pt;
		//calculate the growth direction
		dir.unit();
		vec3d normal = culture->m_pmat->exterior_surface->SurfaceNormal(*surfe, rs[0], rs[1]);
		vec3d gdir = reflect(dir, normal);
		gdir.unit();

		reflSeg.tip(1).pt.r = reflSeg.tip(0).pt.r + (gdir * rdist);

		reflSeg.tip(0).bactive = false;
		reflSeg.tip(1).bactive = true;

		//set seed and vessel to match
		reflSeg.m_nvessel = seg.m_nvessel;
		reflSeg.seed(seg.seed());
		reflSeg.tip(0).nseed = seg.seed();
		reflSeg.tip(1).nseed = seg.seed();
		reflSeg.tip(0).nvessel = seg.m_nvessel;
		reflSeg.tip(1).nvessel = seg.m_nvessel;
		reflSeg.tip(0).connected = &seg;
		reflSeg.Update();
		reflSeg.expected_length = seg.expected_length;

		if (reflSeg.length() >= culture->m_cultParams->min_segment_length)
		{
			//record the index of recent
			culture->AddNewSegmentNoClear(reflSeg);
			found1 = true;
		}
	}
	else
	{

		//probaly too close to the boundary shoot the ray from the other side
		//printf("no segment found\n");
		//the segment is very close to the border of the domain
		//use the direction of the segment to back it off to get a element and normal
		vec3d starting_point = seg.tip(0).pt.r - (segadjdir * 2.0);
		surfe = culture->m_pmat->normal_proj->Project3(starting_point, -dir, rs);
		if (surfe)
		{
			Segment reflSeg;
			reflSeg.tip(0).pt = seg.tip(1).pt;
			//calculate the growth direction
			dir.unit();
			vec3d normal = culture->m_pmat->exterior_surface->SurfaceNormal(*surfe, rs[0], rs[1]);
			vec3d gdir = reflect(dir, normal);
			gdir.unit();

			reflSeg.tip(1).pt.r = reflSeg.tip(0).pt.r + (gdir * rdist);

			reflSeg.tip(0).bactive = false;
			reflSeg.tip(1).bactive = true;

			//set seed and vessel to match
			reflSeg.m_nvessel = seg.m_nvessel;
			reflSeg.seed(seg.seed());
			reflSeg.tip(0).nseed = seg.seed();
			reflSeg.tip(1).nseed = seg.seed();
			reflSeg.tip(0).nvessel = seg.m_nvessel;
			reflSeg.tip(1).nvessel = seg.m_nvessel;
			reflSeg.tip(0).connected = &seg;
			reflSeg.Update();
			reflSeg.expected_length = seg.expected_length;

			if (reflSeg.length() >= culture->m_cultParams->min_segment_length)
			{
				culture->AddNewSegmentNoClear(reflSeg);
				found1 = true;
			}
		}
	}
	auto rseg = culture->RecentSegments();
	if (rseg.size())
	{
		size_t sz = rseg.size();
		rseg[sz - 1]->tip(1).bactive = true;
	}
	//see if the previous tip needs reactivated
}
MBC::MBC(FEModel * model) : FEMaterial(model) //-V730
{
	
}

bool MBC::acceptBoundary(FEAngioMaterialBase* mat0, FEAngioMaterialBase* mat1)
{
	assert(mat0 == culture->m_pmat);
	return (mat0->GetCommonAngioProperties()->bc->angio_boundary_groups & mat1->GetCommonAngioProperties()->bc->angio_boundary_groups) != 0;
}

PassThroughMBC::PassThroughMBC(FEModel* model) : MBC(model)
{
	
}


//this function splits the segment between the cultures mat0 is the originating material while 
//mat1 is the new material
void PassThroughMBC::handleBoundary(FEAngioMaterialBase* mat0, FEAngioMaterialBase* mat1, Segment & seg)
{
	GridPoint oldtip = seg.tip(1).pt;
	vec3d dir = seg.tip(1).pt.r - seg.tip(0).pt.r;
	FESurface * surf = mat0->exterior_surface;
	double tdist = seg.length();
	FESolidElement & se = seg.tip(0).pt.ndomain->Element(seg.tip(0).pt.elemindex);
	std::vector<int> & edinices = mat0->m_pangio->m_fe_element_data[se.GetID()].surfacesIndices;
	mat0->normal_proj->SetSearchRadius(seg.length() * 2);
	double rs[3];
	double g;
	
	bool found0 = false;
	for (size_t i = 0; i < edinices.size(); i++)
	{
		FESurfaceElement & surfe = reinterpret_cast<FESurfaceElement&>(surf->ElementRef(mat0->m_pangio->m_fe_element_data[se.GetID()].surfacesIndices[i]));
		if (mat0->exterior_surface->Intersect(surfe, seg.tip(0).pt.r, -dir, rs, g, 0.0001))//see the epsilon in FIndSolidElement
		{
			//set last_goodpt
			seg.tip(1).pt.r = surf->Local2Global(surfe, rs[0], rs[1]);
			seg.tip(1).pt.q = mat0->m_pangio->FindRST(seg.tip(1).pt.r, vec2d(rs[0], rs[1]), &se);
			seg.tip(1).pt.elemindex = seg.tip(0).pt.elemindex;
			seg.tip(1).pt.nelem = seg.tip(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			seg.tip(1).bactive = false;
			seg.Update();
			
			if (seg.length() > mat0->m_cultureParams.min_segment_length)
			{
				seg.SetFlagOn(Segment::BC_DEAD);
				mat0->m_cult->AddSegment(seg);
				found0 = true;
				break;
			}
			else
			{
				seg.SetFlagOff(Segment::BC_DEAD);
			}
				
		}
	}
	if (!found0)
	{
		FESurfaceElement * surfe = culture->m_pmat->normal_proj->Project(seg.tip(0).pt.r, -dir, rs);
		if (surfe)
		{
			int eindex = surfe->m_elem[0];

			vec3d pos = surf->Local2Global(*surfe, rs[0], rs[1]);
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			seg.tip(1).pt.elemindex = eindex - culture->m_pmat->meshOffsets.at(seg.tip(1).pt.ndomain);
			seg.tip(1).pt.nelem = eindex + 1;

			FESolidElement & se = reinterpret_cast<FESolidElement&>(seg.tip(1).pt.ndomain->ElementRef(seg.tip(1).pt.elemindex));
			seg.tip(1).pt.r = surf->Local2Global(*surfe, rs[0], rs[1]);
			seg.tip(1).pt.q = mat0->m_pangio->FindRST(seg.tip(1).pt.r, vec2d(rs[0], rs[1]), &se);
			seg.tip(1).bactive = false;
			seg.Update();
			
			if ((seg.length() > mat0->m_cultureParams.min_segment_length) && (seg.length() < tdist))
			{
				seg.SetFlagOn(Segment::BC_DEAD);
				mat0->m_cult->AddSegment(seg);
				found0 = true;
			}
			else
			{
				seg.SetFlagOff(Segment::BC_DEAD);
			}

		}
	}
	
	
	//project from the opposite direction
	bool found1 = false;
	Segment s2= seg;//copies the seed and segment number 
	s2.SetFlagOff(Segment::BC_DEAD);
	s2.tip(1).pt = oldtip;
	s2.tip(0) = seg.tip(1);
	dir = -dir;
	
	surf = mat1->exterior_surface;
	
	FESolidElement & se1 = oldtip.ndomain->Element(oldtip.elemindex);
	
	std::vector<int> & edinices1 = mat1->m_pangio->m_fe_element_data[se1.GetID()].surfacesIndices;
	s2.Update();
	mat1->normal_proj->SetSearchRadius(s2.length() * 2);
	for (size_t i = 0; i < edinices1.size(); i++)
	{
		FESurfaceElement & surfe = reinterpret_cast<FESurfaceElement&>(surf->ElementRef(mat1->m_pangio->m_fe_element_data[se1.GetID()].surfacesIndices[i]));
		if (mat1->exterior_surface->Intersect(surfe, oldtip.r, -dir, rs, g, 0.0001))//see the epsilon in FIndSolidElement
		{
			//set last_goodpt
			s2.tip(0).pt.r = surf->Local2Global(surfe, rs[0], rs[1]);
			s2.tip(0).pt.elemindex = oldtip.elemindex;
			s2.tip(0).pt.nelem = oldtip.nelem;
			s2.tip(0).pt.ndomain = oldtip.ndomain;
			s2.tip(0).pt.q = mat1->m_pangio->FindRST(s2.tip(0).pt.r, vec2d(rs[0], rs[1]), &se1);
			
			s2.tip(0).bactive = false;
			s2.tip(1).bactive = true;
			
			s2.tip(0).pt.r = mat1->m_pangio->Position(s2.tip(0).pt);
			s2.Update();
			if (s2.length() > mat1->m_cultureParams.min_segment_length)
			{
				mat1->m_cult->AddNewSegment(s2);
				found1 = true;
				break;
			}
			else
			{
				seg.SetFlagOff(Segment::BC_DEAD);
			}
		}
	}
	if (!found1)
	{
		FESurfaceElement * surfe = mat1->normal_proj->Project(s2.tip(1).pt.r, -dir, rs);
		if (surfe)
		{
			int eindex = surfe->m_elem[0];

			vec3d pos = surf->Local2Global(*surfe, rs[0], rs[1]);
			s2.tip(0).pt.ndomain = s2.tip(1).pt.ndomain;
			s2.tip(0).pt.elemindex = eindex - mat1->meshOffsets.at(s2.tip(1).pt.ndomain);
			s2.tip(0).pt.nelem = eindex + 1;

			FESolidElement & se = reinterpret_cast<FESolidElement&>(s2.tip(1).pt.ndomain->ElementRef(s2.tip(0).pt.elemindex));
			s2.tip(0).pt.r = surf->Local2Global(*surfe, rs[0], rs[1]);
			s2.tip(0).pt.q = mat1->m_pangio->FindRST(s2.tip(0).pt.r, vec2d(rs[0], rs[1]), &se);

			s2.Update();
			if ((s2.length() > mat1->m_cultureParams.min_segment_length) && (s2.length() < tdist))
			{
				s2.tip(1).bactive = true;
				s2.tip(0).bactive = false;
				//mat1's recents may need to be cleared
				mat1->m_cult->AddSegment(s2);
				found1 = true;
			}
			else
			{
				seg.SetFlagOff(Segment::BC_DEAD);
			}

		}
	}

	//add a handler for when the endpoints aren't close
	if ((found1 && found0) || (!found0 && found1))
	{
		if ((mat0->m_cultureParams.min_segment_length > 0.0) && ((seg.tip(1).pt.r - s2.tip(0).pt.r).norm() >(2.0 * mat0->m_cultureParams.min_segment_length)))
		{
			printf("warning growing entirely through a material\n");
			mat0->m_pangio->m_fe_element_data[se.GetID()].flags |= 1;
			
		}
	}
	else if ((found0 && !found1))
	{
		printf("only one segment found\n");
		mat0->m_pangio->m_fe_element_data[se.GetID()].flags |= 2;
	}
}
