// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

#include "FECore/FECoreKernel.h"
#include "AngioFETask.h"
#include "FECore/FECoreFactory.h"
#include "FESproutBodyForce.h"
#include "FEAngioMaterial.h"
#include "AngioPlot.h"
#ifdef SVN
#include "svnrev.h"
#else
#define SVNREVISION 0
#endif

//-----------------------------------------------------------------------------
FEPluginFactory_T<AngioFETask                , FETASK_ID    > angiofe_task_factory         ("angio"                 );
FEPluginFactory_T<FESproutBodyForce          , FEBODYLOAD_ID> angio_sprout_factory         ("sprout"                );
FEPluginFactory_T<FEAngioMaterial            , FEMATERIAL_ID> angio_mat_factory            ("angio"                 );
FEPluginFactory_T<FEPressureMaterial         , FEMATERIAL_ID> pressure_mat_factory         ("pressure"              );
FEPluginFactory_T<FEPlotAngioStress          , FEPLOTDATA_ID> plot_angio_stress            ("angio stress"          );
FEPluginFactory_T<FEPlotVesselStress         , FEPLOTDATA_ID> plot_vessel_stress           ("vessel stress"         );
FEPluginFactory_T<FEPlotMatrixStress         , FEPLOTDATA_ID> plot_matrix_stress           ("matrix stress"         );
FEPluginFactory_T<FEPlotVesselWeight         , FEPLOTDATA_ID> plot_vessel_weight           ("vessel weight"         );
FEPluginFactory_T<FEPlotMatrixWeight         , FEPLOTDATA_ID> plot_matrix_weight           ("matrix weight"         );
FEPluginFactory_T<FEPlotAngioEffectiveStress , FEPLOTDATA_ID> plot_angio_eff_stress        ("angio effective stress");
FEPluginFactory_T<FEPlotAngioCollagenFibers  , FEPLOTDATA_ID> plot_angio_fibers            ("angio collagen fiber"  );
FEPluginFactory_T<FEPlotAngioECMDensity      , FEPLOTDATA_ID> plot_angio_ecm               ("angio ECM density"     );
FEPluginFactory_T<FEPlotAngioECMAlpha        , FEPLOTDATA_ID> plot_angio_alpha             ("angio ECM alpha"       );
FEPluginFactory_T<FEPlotAngioGradient        , FEPLOTDATA_ID> plot_angio_gradient          ("angio gradient"        );
FEPluginFactory_T<FEPlotAngioGradientCenter  , FEPLOTDATA_ID> plot_angio_gradient_center   ("angio gradient center" );
FEPluginFactory_T<FEPlotAngioMaterialHop     , FEPLOTDATA_ID> plot_angio_material_hop      ("angio material hop"    );
FEPluginFactory_T<FEPlotAngioSegmentBadGrowth, FEPLOTDATA_ID> plot_angio_segment_bad_growth("angio bad growth"      );

//-----------------------------------------------------------------------------
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
}

FECORE_EXPORT void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 1;
	minor = 0;
	patch = SVNREVISION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT int PluginNumClasses()
{
	return 17;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT FECoreFactory * PluginGetFactory(int i)
{
	switch (i)
	{
	case 0:
		return &angiofe_task_factory;
	case 1:
		return &angio_sprout_factory;
	case 2:
		return &angio_mat_factory;
	case 3:
		return &pressure_mat_factory;
	case 4:
		return &plot_angio_stress;
	case 5:
		return &plot_angio_eff_stress;
	case 6:
		return &plot_angio_fibers;
	case 7:
		return &plot_angio_ecm;
	case 8:
		return &plot_angio_alpha;
	case 9:
		return &plot_angio_gradient;
	case 10:
		return &plot_angio_gradient_center;
	case 11:
		return &plot_angio_material_hop;
	case 12:
		return &plot_angio_segment_bad_growth;
	case 13:
		return &plot_vessel_stress;
	case 14:
		return &plot_matrix_stress;
	case 15:
		return &plot_vessel_weight;
	case 16:
		return &plot_matrix_weight;
	default:
		return nullptr;
	}
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginCleanup()
{

}
