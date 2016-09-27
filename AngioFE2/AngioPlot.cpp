#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"

extern FEAngio* pfeangio;

//-----------------------------------------------------------------------------
bool FEPlotAngioStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->FindAngioMaterial(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			mat3ds sj = pmat->AngioStress(*(FEAngioMaterialPoint::FindAngioMaterialPoint(&mp)));
			s += sj;
		}
		s /= (double) nint;

		str << s;
	}
	return true;
};

//-----------------------------------------------------------------------------
bool FEPlotAngioEffectiveStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = dynamic_cast<FEAngioMaterial*>(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			mat3ds sj = pmat->Stress(*(FEAngioMaterialPoint::FindAngioMaterialPoint(&mp)));
			s += pt.m_s - sj;
		}
		s /= (double) nint;

		str << s;
	}
	return true;
};

//-----------------------------------------------------------------------------
bool FEPlotAngioCollagenFibers::Save(FEMesh& m, FEDataStream& a)
{
	if (pfeangio == 0) return false;
	//multiple materials average their *
	for (int i = 0; i < pfeangio->m_fem.GetMesh().Nodes(); i++)
	{
		if (pfeangio->m_fe_node_data.count(pfeangio->m_fem.GetMesh().Node(i).GetID()))
		{
			a << pfeangio->m_fe_node_data[pfeangio->m_fem.GetMesh().Node(i).GetID()].m_collfib;
		}
		else
		{
			a << vec3d(0, 0, 0);//is all zero's okay for this parameter
		}
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotAngioECMDensity::Save(FEMesh& m, FEDataStream& a)
{

	//multiple materials average their ecm denisties so this reflects that
	for (int i = 0; i < pfeangio->m_fem.GetMesh().Nodes(); i++)
	{
		if (pfeangio->m_fe_node_data.count(pfeangio->m_fem.GetMesh().Node(i).GetID()))
		{
			a << pfeangio->m_fe_node_data[pfeangio->m_fem.GetMesh().Node(i).GetID()].m_ecm_den;
		}
		else
		{
			a << 0.0;//is zero okay for this parameter
		}
	}

	return true;
}

bool FEPlotAngioECMAlpha::Save(FEMesh& m, FEDataStream& a)
{
	if (pfeangio == 0) return false;
	//multiple materials average their *
	for (int i = 0; i < pfeangio->m_fem.GetMesh().Nodes(); i++)
	{
		if (pfeangio->m_fe_node_data.count(pfeangio->m_fem.GetMesh().Node(i).GetID()))
		{
			a << pfeangio->m_fe_node_data[pfeangio->m_fem.GetMesh().Node(i).GetID()].alpha;
		}
		else
		{
			a << 0.0;//is zero okay for this parameter
		}
	}

	return true;
}