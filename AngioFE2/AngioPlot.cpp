#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include <FECore/FESolidDomain.h>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"
#include "FECore/FEMesh.h"
#include <unordered_map>
#include <algorithm>
#include <FEBioMech/FEViscoElasticMaterial.h>
#include <FEBioMix/FEMultiphasic.h>
#include "FEBioMix/FEBiphasicSolute.h"
#include "FEBioMix/FETriphasic.h"
#include "FEBioMix/FEMultiphasicSolidDomain.h"
#include "FEBioMix/FEMultiphasicShellDomain.h"

extern FEAngio* pfeangio;

//-----------------------------------------------------------------------------
bool FEPlotAngioStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
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
			FEAngioMaterialPoint * angio_mp = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			mat3ds & sj = angio_mp->m_as;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEElasticMaterialPoint& matrix_elastic = *angioPt->matPt->ExtractData<FEElasticMaterialPoint>();
			mat3ds & sj = matrix_elastic.m_s;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotVesselStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEElasticMaterialPoint& vessel_elastic = *angioPt->vessPt->ExtractData<FEElasticMaterialPoint>();
			mat3ds & sj = vessel_elastic.m_s;
			
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotVesselWeight::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		double s = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			s += angioPt->vessel_weight;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixWeight::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		double s = 0.0;
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			s += angioPt->matrix_weight;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixTangent::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		tens4ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			tens4ds ten = pmat->GetMatrixMaterial()->GetElasticMaterial()->Tangent(mp);
			tens4ds sj = ten;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixViscoStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEViscoElasticMaterialPoint& matrix_visco_elastic = *angioPt->matPt->ExtractData<FEViscoElasticMaterialPoint>();
			
			mat3ds sj =  matrix_visco_elastic.m_se;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixElasticStress::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint& mp = *(el.GetMaterialPoint(j));
			FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp);
			FEElasticMaterialPoint& emp = *angioPt->matPt->Next()->ExtractData<FEElasticMaterialPoint>();

			mat3ds sj = emp.m_s;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixElastic_m_Q::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3d s;
		s.zero();
		for (int j = 0; j<nint; ++j)
		{
			FEMaterialPoint * mp = (el.GetMaterialPoint(j));
			FEElasticMaterialPoint*  emp = mp->ExtractData<FEElasticMaterialPoint>();

			mat3d sj = emp->m_Q;
			s += sj;
		}
		s /= static_cast<double>(nint);

		str << s;
	}
	return true;
}


FEPlotAngioGradientCenter::FEPlotAngioGradientCenter(FEModel* pfem): FEDomainData(PLT_VEC3F, FMT_ITEM)
{
};

bool FEPlotAngioGradientCenter::Save(FEDomain& d, FEDataStream& str)
{
	//note the interpolation in postview may make this look incorrect
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
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
FEPlotAngioMaterialHop::FEPlotAngioMaterialHop(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
};

bool FEPlotAngioMaterialHop::Save(FEDomain& d, FEDataStream& str)
{
	//note the interpolation in postview may make this look incorrect
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);

		float val =  static_cast<float>((pmat->m_pangio->m_fe_element_data[el.GetID()].flags & 1));

		str << val;
	}
	return true;
};
FEPlotAngioSegmentBadGrowth::FEPlotAngioSegmentBadGrowth(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM)
{
};

bool FEPlotAngioSegmentBadGrowth::Save(FEDomain& d, FEDataStream& str)
{
	//note the interpolation in postview may make this look incorrect
	FEAngioMaterial* pmat = pfeangio->GetAngioComponent(d.GetMaterial());
	if (pmat == nullptr) return false;

	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i = 0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);

		float val = static_cast<float>((pmat->m_pangio->m_fe_element_data[el.GetID()].flags & 2));

		str << val;
	}
	return true;
};


bool FEPlotAngioGradient::Save(FEMesh & m, FEDataStream & a)
{
	//this has problems on the boundaries between materials
	std::unordered_map<int, vec3d> gradients;
	if (pfeangio == nullptr) return false;
	//
	FEMesh & mesh = pfeangio->m_fem->GetMesh();
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


		
		for (size_t i = 0; i < se.m_node.size(); i++)
		{
			vec3d pt = vec3d(nr[i], ns[i], nt[i]);
			FENode & node = mesh.Node(se.m_node[i]);
			
			if (pfeangio->GetAngioComponent(d.GetMaterial())->Overwrite() || (gradients.count(node.GetID()) == 0 ))
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
	for (int i = 0; i < pfeangio->m_fem->GetMesh().Nodes(); i++)
	{
		if (pfeangio->m_fe_node_data.count(pfeangio->m_fem->GetMesh().Node(i).GetID()))
		{
			a << pfeangio->m_fe_node_data[pfeangio->m_fem->GetMesh().Node(i).GetID()].m_ecm_den;
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
	for (int i = 0; i < pfeangio->m_fem->GetMesh().Nodes(); i++)
	{
		if (pfeangio->m_fe_node_data.count(pfeangio->m_fem->GetMesh().Node(i).GetID()))
		{
			a << pfeangio->m_fe_node_data[pfeangio->m_fem->GetMesh().Node(i).GetID()].alpha;
		}
		else
		{
			a << 0.0;//is zero okay for this parameter
		}
	}

	return true;
}

//stopgap plot implementation

//-----------------------------------------------------------------------------
// find the local solute ID, given a global ID. If the material is not a 
// biphasic-solute, triphasic, or multiphasic material, this returns -1.
int GetLocalSoluteID(FEMaterial* pm, int nsol)
{
	// figure out the solute ID to export. This depends on the material type.
	int nsid = -1;
	FEBiphasicSolute* psm = dynamic_cast<FEBiphasicSolute*> (pm);
	if (psm)
	{
		// Check if this solute is present in this specific biphasic-solute mixture
		bool present = (psm->GetSolute()->GetSoluteID() == nsol);
		if (!present) return false;
		nsid = 0;
	}

	FETriphasic* ptm = dynamic_cast<FETriphasic*> (pm);
	if (ptm)
	{
		// Check if this solute is present in this specific triphasic mixture
		if (ptm->m_pSolute[0]->GetSoluteID() == nsol) nsid = 0;
		else if (ptm->m_pSolute[1]->GetSoluteID() == nsol) nsid = 1;
	}

	FEMultiphasic* pmm = dynamic_cast<FEMultiphasic*> (pm);
	if (pmm)
	{
		// Check if this solute is present in this specific multiphasic mixture
		for (int i = 0; i<pmm->Solutes(); ++i)
			if (pmm->GetSolute(i)->GetSoluteID() == nsol) { nsid = i; break; }
	}
	return nsid;
}
bool FEPlotMatrixConectrationGradient::Save(FEDomain& d, FEDataStream& str)
{
	FEAngioMaterial * angio_mat = dynamic_cast<FEAngioMaterial*>(d.GetMaterial());
	int nsid = GetLocalSoluteID(angio_mat->GetMatrixMaterial(), 1);

	// make sure we have a valid index
	if (nsid == -1) return false;

	int N = d.Elements();
	for (int i = 0; i<N; ++i)
	{
		FEElement& el = d.ElementRef(i);

		// calculate average concentration
		double ew = 0;
		for (int j = 0; j<el.GaussPoints(); ++j)
		{
			FEMaterialPoint& mp = *el.GetMaterialPoint(j);
			FESolutesMaterialPoint* pt = (mp.ExtractData<FESolutesMaterialPoint>());

			if (pt) ew += pt->m_ca[nsid];
		}

		ew /= el.GaussPoints();

		str << ew;
	}
	return true;
}

//-----------------------------------------------------------------------------
bool FEPlotMatrixSBMConectration::Save(FEDomain &d, FEDataStream& str)
{
	int i, j;
	double ew;
	const int sbm_id = 1;
	FEMultiphasicSolidDomain* pmd = dynamic_cast<FEMultiphasicSolidDomain*>(&d);
	FEMultiphasicShellDomain* psd = dynamic_cast<FEMultiphasicShellDomain*>(&d);
	FEAngioMaterial * angio_mat = dynamic_cast<FEAngioMaterial*>(d.GetMaterial());
	assert(angio_mat);
	FEMaterial * matrix_mat = angio_mat->GetMatrixMaterial();
	if (pmd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*>(matrix_mat);
		assert(pm);
		// Check if this solid-bound molecule is present in this specific multiphasic mixture
		int sid = -1;
		for (i = 0; i<pm->SBMs(); ++i)
			if (pm->GetSBM(i)->GetSBMID() == sbm_id) { sid = i; break; }
		if (sid == -1) return false;

		for (i = 0; i<pmd->Elements(); ++i)
		{
			FESolidElement& el = pmd->Element(i);

			// calculate average concentration
			ew = 0;
			for (j = 0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());

				if (st) ew += pm->SBMConcentration(mp, sid);
			}

			ew /= el.GaussPoints();

			str << ew;
		}
		return true;
	}
	else if (psd)
	{
		FEMultiphasic* pm = dynamic_cast<FEMultiphasic*> (d.GetMaterial());
		// Check if this solid-bound molecule is present in this specific multiphasic mixture
		int sid = -1;
		for (i = 0; i<pm->SBMs(); ++i)
			if (pm->GetSBM(i)->GetSBMID() == sbm_id) { sid = i; break; }
		if (sid == -1) return false;

		for (i = 0; i<psd->Elements(); ++i)
		{
			FEShellElement& el = psd->Element(i);

			// calculate average concentration
			ew = 0;
			for (j = 0; j<el.GaussPoints(); ++j)
			{
				FEMaterialPoint& mp = *el.GetMaterialPoint(j);
				FESolutesMaterialPoint* st = (mp.ExtractData<FESolutesMaterialPoint>());

				if (st) ew += pm->SBMConcentration(mp, sid);
			}

			ew /= el.GaussPoints();

			str << ew;
		}
		return true;
	}
	return false;
}