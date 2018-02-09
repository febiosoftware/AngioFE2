#include "StdAfx.h"
#include "ECMInitializer.h"
#include <FECore/FEMesh.h>
#include "FEAngioMaterialPoint.h"
#include "FEAngioMaterialBase.h"
#include "FEAngio.h"

void ECMInitializerConstant::seedECMDensity(FEAngioMaterialBase * mat)
{
	std::vector<int> matls;
	matls.emplace_back(mat->GetID_ang());

	mat->m_pangio->ForEachElement([&](FESolidElement & el, FESolidDomain & dom)
	{
		for (int i = 0; i< el.GaussPoints();i++)
		{
			FEAngioMaterialPoint *mp = FEAngioMaterialPoint::FindAngioMaterialPoint(el.GetMaterialPoint(i));
			assert(mp);
			mp->ref_ecm_density = mat->m_cultureParams.m_matrix_density;
			mp->anisotropy = mat->GetAnisotropy();
		}
	}, matls);
}