#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include <FECore/FESolidDomain.h>
#include "FEAngioMaterial.h"
#include "Grid.h"

extern FEAngio* pfeangio;

//-----------------------------------------------------------------------------
bool FEPlotAngioStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = dynamic_cast<FEAngioMaterial*>(d.GetMaterial());
	if (pmat == 0) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s(0.0);
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			mat3ds sj = pmat->Stress(mp);
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
	if (pmat == 0) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s(0.0);
		for (int j=0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEElasticMaterialPoint& pt = *(mp.ExtractData<FEElasticMaterialPoint>());
			mat3ds sj = pmat->Stress(mp);
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

	Grid& grid = pfeangio->GetGrid();
	assert(grid.Nodes() == m.Nodes());

	int NN = grid.Nodes();
	for (int i=0; i<NN; ++i) a << grid.GetNode(i).m_collfib;

	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotAngioECMDensity::Save(FEMesh& m, FEDataStream& a)
{
	if (pfeangio == 0) return false;

	Grid& grid = pfeangio->GetGrid();
	assert(grid.Nodes() == m.Nodes());

	int NN = grid.Nodes();
	for (int i=0; i<NN; ++i) a << grid.GetNode(i).m_ecm_den;

	return true;
}
