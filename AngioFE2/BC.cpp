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
BC::BC(FEAngio& angio) : m_angio(angio)
{
}                                   

//-----------------------------------------------------------------------------
BC::~BC()
{

}

//-----------------------------------------------------------------------------
// Checks if a newly-created segment violates the boundary faces of the element in which it occupies
void BC::CheckBC(Segment &seg, Culture * culture)
{
	
	//new implementation may be run on all segments will add the segment once the boundaries are safe
	//remove the check of whether or not to check the boundary condition in Culture
	assert(seg.tip(0).pt.nelem != -1);
	//make sure that the defining tip has not drifted too far
	vec3d pq = m_angio.Position(seg.tip(0).pt);
	assert((m_angio.Position(seg.tip(0).pt) - seg.tip(0).pt.r).norm() < 1.0);
	//if the second segment is unititialized make sure it is not in the same element as tip(0)
	FESolidElement * tse = dynamic_cast<FESolidElement*>(&seg.tip(0).pt.ndomain->ElementRef(seg.tip(0).pt.elemindex));

	if (seg.tip(1).pt.nelem == -1)
	{
		
		double arr[3];
		if (tse && m_angio.IsInsideHex8(tse,seg.tip(1).pt.r, m_angio.GetMesh() , arr))
		{
			//just copy the data from tip(0)
			seg.tip(1).pt.q.x = arr[0]; seg.tip(1).pt.q.y = arr[1]; seg.tip(1).pt.q.z = arr[2];
			seg.tip(1).pt.nelem = seg.tip(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			seg.tip(1).pt.elemindex = seg.tip(0).pt.elemindex;
		}
	}
	//this is the costliest part of the boundary check if reached
	if (seg.tip(1).pt.nelem == -1)
	{
		culture->m_pmat->FindGridPoint(seg.tip(1).pt.r ,seg.tip(1).pt);
	}
	//check if the segment is in another angio material
	if (seg.tip(1).pt.nelem == -1)
	{
		FEMesh * mesh = m_angio.GetMesh();
		double r[3];
		FESolidElement * se = mesh->FindSolidElement(seg.tip(1).pt.r, r);
		FEAngioMaterial * angm;
		if (se && (angm = dynamic_cast<FEAngioMaterial*>(m_angio.m_fem.GetMaterial(se->GetMatID()))))
		{
			printf("growing into angio material\n");
			seg.SetFlagOn(Segment::BC_DEAD);
			return;
		}
	}
	
	//if both are in the same element just add the segment
	if (seg.tip(0).pt.nelem == seg.tip(1).pt.nelem)
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
	if (seg.tip(1).pt.ndomain != nullptr && seg.tip(1).pt.elemindex >= 0)
	{
		FEElement & t1se = seg.tip(1).pt.ndomain->ElementRef(seg.tip(1).pt.elemindex);
		FEElement & t0se = seg.tip(0).pt.ndomain->ElementRef(seg.tip(0).pt.elemindex);
		if (t1se.GetMatID() == t0se.GetMatID())
		{
			if (seg.length() >= culture->m_cultParams->min_segment_length)
			{
				culture->AddSegment(seg);
				return;
			}
		}
	}
	
	vec3d dir = seg.tip(1).pt.r - seg.tip(0).pt.r;
	double dist = dir.norm();
	assert(dir.norm() > 0.0);

	dir.unit();

	FESolidDomain & sd = reinterpret_cast<FESolidDomain&>(*seg.tip(0).pt.ndomain);
	FESolidElement * se = dynamic_cast<FESolidElement*>(&sd.ElementRef(seg.tip(0).pt.elemindex));
	assert(se);//make sure we have an element

	FEMaterial * mat = sd.GetMaterial();
	vec3d lastgood_pt;
	double rs[2];
	double g;
	assert(culture->m_pmat);
	FESurface * surf = culture->m_pmat->exterior_surface;
	assert(m_angio.m_fe_element_data[se->GetID()].surfacesIndices.size());
	int hcount = 0;//counts how many times the boundary is handled
	std::vector<int> & edinices = m_angio.m_fe_element_data[se->GetID()].surfacesIndices;
	culture->m_pmat->normal_proj->SetSearchRadius(seg.length() * 2);
	for (size_t i = 0; i < edinices.size(); i++)
	{
		FESurfaceElement & surfe = reinterpret_cast<FESurfaceElement&>(surf->ElementRef(m_angio.m_fe_element_data[se->GetID()].surfacesIndices[i]));
		if (culture->m_pmat->exterior_surface->Intersect(surfe, seg.tip(0).pt.r, -dir, rs, g, 0.0001))//see the epsilon in FIndSolidElement
		{
			//set last_goodpt
			lastgood_pt = surf->Local2Global(surfe, rs[0], rs[1]);

			seg.tip(1).pt.elemindex = seg.tip(0).pt.elemindex;
			seg.tip(1).pt.nelem = seg.tip(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			return HandleBoundary(seg, culture, lastgood_pt, rs, se);
		}
		
	}
	FESurfaceElement * surfe = culture->m_pmat->normal_proj->Project(seg.tip(0).pt.r, -dir, rs);
	if (!surfe)
	{
		//printf("no surface element found\n");
		//assert(false);
		surfe = culture->m_pmat->normal_proj->Project3(seg.tip(0).pt.r + (dir * 2.0), -dir, rs);
		if (surfe)
		{
			lastgood_pt = surf->Local2Global(*surfe, rs[0], rs[1]);

			seg.tip(1).pt.elemindex = seg.tip(0).pt.elemindex;
			seg.tip(1).pt.nelem = seg.tip(0).pt.nelem;
			seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
			return HandleBoundary(seg, culture, lastgood_pt, rs, se);
		}
		else
		{
			printf("segment adjustment failed\n");
		}
	}
	else
	{
		int eindex = surfe->m_elem[0];
			
		vec3d pos = surf->Local2Global(*surfe, rs[0], rs[1]);
		seg.tip(1).pt.ndomain = seg.tip(0).pt.ndomain;
		seg.tip(1).pt.elemindex = eindex - culture->m_pmat->meshOffsets[seg.tip(1).pt.ndomain];
		seg.tip(1).pt.nelem = eindex + 1;
		
		FESolidElement & se = reinterpret_cast<FESolidElement&>(seg.tip(1).pt.ndomain->ElementRef(seg.tip(1).pt.elemindex));
		return HandleBoundary(seg, culture, pos, rs, &se);
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
	FESolidElement * se = dynamic_cast<FESolidElement*>(&sd.ElementRef(seg.tip(0).pt.elemindex));
	assert(se);//make sure we have an element


	return false;

}
void StopBC::HandleBoundary(Segment & seg, Culture * culture, vec3d lastGoodPt, double * rs, FESolidElement * se)
{
	//fill in the pt's data and add the segment
	//remaining distance is ignored
	FEMesh * mesh = m_angio.GetMesh();
	seg.tip(1).pt.r = lastGoodPt;
	seg.SetFlagOn(Segment::BC_DEAD);
	if (culture->m_pmat->FindGridPoint(lastGoodPt, seg.tip(1).pt.ndomain, seg.tip(1).pt.elemindex, seg.tip(1).pt))
	{
		assert(seg.tip(1).pt.nelem != -1);
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
		else
		{
			
		}
	}

	else
	{
		//this denotes a mismatch between the FindSolidElement and FESurface output
		//printf("boundary segment out of bounds attempting to correct this\n");
		vec2d nrs(rs[0], rs[1]);
		seg.tip(1).pt.q = m_angio.FindRST(lastGoodPt, nrs, se);
		seg.tip(1).pt.r = lastGoodPt;
		seg.tip(1).pt.nelem = se->GetID();
		seg.tip(1).pt.ndomain = culture->m_pmat->domainptrs[0]; //hack
		seg.tip(1).pt.elemindex = se->GetID() - 1 - culture->m_pmat->meshOffsets[seg.tip(1).pt.ndomain];
		seg.tip(1).pt.r = m_angio.Position(seg.tip(1).pt);
		seg.Update();
		seg.tip(0).bactive = false;
		seg.tip(1).bactive = false;
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			culture->AddSegment(seg);
		}
		else
		{
			printf("boundary segment dropped due to zero length\n");
			//assert(false);
		}
		
	}
}
void BouncyBC::HandleBoundary(Segment & seg, Culture * culture, vec3d lastGoodPt, double * rs, FESolidElement * se)
{
	//fill in the pt's data and add the segment
	//remaining distance is ignored
	double rdist = (lastGoodPt - seg.tip(1).pt.r).norm();
	FEMesh * mesh = m_angio.GetMesh();
	assert(rdist >= 0.0);
	vec3d prevhead = seg.tip(1).pt.r;
	vec3d segadjdir = seg.tip(1).pt.r - seg.tip(0).pt.r;
	segadjdir.unit();
	seg.tip(1).pt.r = lastGoodPt;
	seg.SetFlagOn(Segment::BC_DEAD);
	if (culture->m_pmat->FindGridPoint(lastGoodPt, seg.tip(1).pt.ndomain, seg.tip(1).pt.elemindex, seg.tip(1).pt))
	{
		assert(seg.tip(1).pt.nelem != -1);
		seg.Update();
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			seg.tip(0).bactive = false;
			seg.tip(1).bactive = false;
			seg.tip(1).bdyf_id = 1;
			seg.tip(0).bdyf_id = 1;
			culture->AddSegment(seg);
		}
		else
		{
			return;
		}
	}

	else
	{
		//this denotes a mismatch between the FindSolidElement and FESurface output
		//printf("boundary segment out of bounds attempting to correct this\n");
		vec2d nrs(rs[0], rs[1]);
		seg.tip(1).pt.q = m_angio.FindRST(lastGoodPt, nrs, se);
		seg.tip(1).pt.r = lastGoodPt;
		seg.tip(1).pt.nelem = se->GetID();
		seg.tip(1).pt.ndomain = culture->m_pmat->domainptrs[0]; //hack
		seg.tip(1).pt.elemindex = se->GetID() - 1 - culture->m_pmat->meshOffsets[seg.tip(1).pt.ndomain];
		seg.tip(1).pt.r = m_angio.Position(seg.tip(1).pt);
		seg.Update();
		seg.tip(0).bactive = false;
		seg.tip(1).bactive = false;
		seg.tip(1).bdyf_id = 1;
		seg.tip(0).bdyf_id = 1;
		if (seg.length() >= culture->m_cultParams->min_segment_length)
		{
			culture->AddSegment(seg);
		}
		else
		{
			printf("boundary segment dropped due to zero length\n");
			return;
			//assert(false);
		}

	}
	assert(seg.tip(1).pt.nelem != -1);
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
		reflSeg.tip(1).bdyf_id = seg.tip(1).bdyf_id;
		reflSeg.Update();
		reflSeg.tip(1).bdyf_id = 1;
		reflSeg.tip(0).bdyf_id = 1;

		if (reflSeg.length() >= culture->m_cultParams->min_segment_length)
			return culture->AddNewSegment(reflSeg);
		
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
			reflSeg.tip(1).bdyf_id = seg.tip(1).bdyf_id;
			reflSeg.Update();
			reflSeg.tip(1).bdyf_id = 1;
			reflSeg.tip(0).bdyf_id = 1;

			if (reflSeg.length() >= culture->m_cultParams->min_segment_length)
				return culture->AddNewSegment(reflSeg);
		}
		else
		{
			printf("adjustment not working\n");
		}
	}
	
}
