#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"
#include <unordered_map>
#include <algorithm>

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
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
};

bool FEPlotAngioGradientCenter::Save(FEDomain& d, FEDataStream& str)
{
	//note the interpolation in postview may make this look incorrect
	FEAngioMaterial* pmat = pfeangio->FindAngioMaterial(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		
		vec3d zero;
		std::vector<double> densities;
		densities = pfeangio->createVectorOfMaterialParameters(&el, &FEAngioNodeData::m_ecm_den);
		vec3d grad = pfeangio->gradient(&el, densities, zero);

		str << grad;
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
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
};

//-----------------------------------------------------------------------------
bool FEPlotAngioCollagenFibers::Save(FEMesh& m, FEDataStream& a)
{
	if (pfeangio == nullptr) return false;
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
bool FEPlotAngioGradient::Save(FEMesh & m, FEDataStream & a)
{
	//this has problems on the boundaries between materials
	std::unordered_map<int, vec3d> gradients;
	if (pfeangio == nullptr) return false;
	//
	FEMesh & mesh = pfeangio->m_fem.GetMesh();
	pfeangio->ForEachElement([&mesh, &gradients](FESolidElement & se, FESolidDomain & d)
	{
		//these will hold the natural coordinates once the project to nodes is complete 
		double nr[FEElement::MAX_NODES];
		double ns[FEElement::MAX_NODES];
		double nt[FEElement::MAX_NODES];
		//these hold the natural coordinates of the integration points (r,s,t)
		double gr[FEElement::MAX_NODES];
		double gs[FEElement::MAX_NODES];
		double gt[FEElement::MAX_NODES];

			
		for (int i = 0; i < se.Nodes(); i++)
		{
			gr[i] = se.gr(i);
			gs[i] = se.gs(i);
			gt[i] = se.gt(i);
		}

		se.project_to_nodes(gr, nr);
		se.project_to_nodes(gs, ns);
		se.project_to_nodes(gt, nt);
		std::vector<double> densities;
		densities = pfeangio->createVectorOfMaterialParameters(&se, &FEAngioNodeData::m_ecm_den);


		
		for (int i = 0; i < se.m_node.size(); i++)
		{
			vec3d pt = vec3d(nr[i], ns[i], nt[i]);
			FENode & node = mesh.Node(se.m_node[i]);
			
			if (dynamic_cast<FEAngioMaterial*>(d.GetMaterial())->Overwrite() || (gradients.count(node.GetID()) == 0 ))
			{
				gradients[node.GetID()] = pfeangio->gradient(&se, densities, pt);
			}
		}
		
	});


	for (int i = 0; i < mesh.Nodes(); i++)
	{
		if (gradients.count(mesh.Node(i).GetID()))
		{
			a << gradients[mesh.Node(i).GetID()];
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