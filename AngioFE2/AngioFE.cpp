// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

#define FECORE_DLL
#define FECORE_API

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

FEPluginFactory_T<NoFragmentBranching        , FEMATERIAL_ID> no_fragment_branching_factory("no_branch"             );
FEPluginFactory_T<PsuedoDeferedFragmentBranching, FEMATERIAL_ID> psuedo_defered_fragment_branching_factory("psuedo_defered_branch");

FEPluginFactory_T<FENormalDistribution, FEMATERIAL_ID> normal_distribution_factory("normal_distribution");
FEPluginFactory_T<FEExponentialDistribution, FEMATERIAL_ID> exponential_distribution_factory("exponential_distribution");
FEPluginFactory_T<FECauchyDistribution, FEMATERIAL_ID> cauchy_distribution_factory("cauchy_distribution");
FEPluginFactory_T<FEChiSquaredDistribution, FEMATERIAL_ID> chi_squared_distribution_factory("chi_squared_distribution");

FEPluginFactory_T<GrowDirectionModifiers, FEMATERIAL_ID> grow_direction_modifiers_factory("grow_direction_modifiers");

FEPluginFactory_T<DefaultGrowDirectionModifier, FEMATERIAL_ID> default_grow_direction_modifier_factory("default_grow_direction");
FEPluginFactory_T<BranchGrowDirectionModifier, FEMATERIAL_ID> branch_grow_direction_modifier_factory("branch_grow_direction");
FEPluginFactory_T<GradientGrowDirectionModifier, FEMATERIAL_ID> gradient_grow_direction_modifier_factory("gradient_grow_direction");
FEPluginFactory_T<AnastamosisGrowDirectionModifier, FEMATERIAL_ID> anastamosis_grow_direction_modifier_factory("anastamosis_grow_direction");

FEPluginFactory_T<ClassicFragmentSeeder, FEMATERIAL_ID> classic_fragment_seeder_factory("classic");
FEPluginFactory_T<MDByVolumeFragmentSeeder, FEMATERIAL_ID> mdbyvol_fragment_seeder_factory("MDbyVolume");
FEPluginFactory_T<MultiDomainFragmentSeeder, FEMATERIAL_ID> md_fragment_seeder_factory("MD");
FEPluginFactory_T<MDAngVessFileFragmentSeeder, FEMATERIAL_ID> md_file_seeder_factory("from_file");

FEPluginFactory_T<BouncyBC, FEMATERIAL_ID> bouncybc_factory("bouncy");
FEPluginFactory_T<StopBC, FEMATERIAL_ID> stopbc_factory("stop");

FEPluginFactory_T<SameMBC, FEMATERIAL_ID> same_mbc_factory("same");
FEPluginFactory_T<PassThroughMBC, FEMATERIAL_ID> passthrough_mbc_factory("pass_through");

FEPluginFactory_T<FEPlotAngioStress          , FEPLOTDATA_ID> plot_angio_stress            ("angio stress"          );
FEPluginFactory_T<FEPlotVesselStress         , FEPLOTDATA_ID> plot_vessel_stress           ("vessel stress"         );
FEPluginFactory_T<FEPlotMatrixStress         , FEPLOTDATA_ID> plot_matrix_stress           ("matrix stress"         );
FEPluginFactory_T<FEPlotVesselWeight         , FEPLOTDATA_ID> plot_vessel_weight           ("vessel weight"         );
FEPluginFactory_T<FEPlotMatrixWeight         , FEPLOTDATA_ID> plot_matrix_weight           ("matrix weight"         );
FEPluginFactory_T<FEPlotMatrixTangent        , FEPLOTDATA_ID> plot_matrix_tangent          ("matrix tangent"        );
FEPluginFactory_T<FEPlotMatrixViscoStress    , FEPLOTDATA_ID> plot_matrix_visco_stress     ("matrix visco stress"   );
FEPluginFactory_T<FEPlotMatrixElasticStress  , FEPLOTDATA_ID> plot_matrix_elastic_stress   ("matrix elastic stress" );
FEPluginFactory_T<FEPlotAngioEffectiveStress , FEPLOTDATA_ID> plot_angio_eff_stress        ("angio effective stress");
FEPluginFactory_T<FEPlotAngioCollagenFibers  , FEPLOTDATA_ID> plot_angio_fibers            ("angio collagen fiber"  );
FEPluginFactory_T<FEPlotAngioECMDensity      , FEPLOTDATA_ID> plot_angio_ecm               ("angio ECM density"     );
FEPluginFactory_T<FEPlotAngioECMAlpha        , FEPLOTDATA_ID> plot_angio_alpha             ("angio ECM alpha"       );
FEPluginFactory_T<FEPlotAngioGradient        , FEPLOTDATA_ID> plot_angio_gradient          ("angio gradient"        );
FEPluginFactory_T<FEPlotAngioGradientCenter  , FEPLOTDATA_ID> plot_angio_gradient_center   ("angio gradient center" );
FEPluginFactory_T<FEPlotAngioMaterialHop     , FEPLOTDATA_ID> plot_angio_material_hop      ("angio material hop"    );
FEPluginFactory_T<FEPlotAngioSegmentBadGrowth, FEPLOTDATA_ID> plot_angio_segment_bad_growth("angio bad growth"      );

extern "C"
{

//-----------------------------------------------------------------------------
FECORE_EXPORT  unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
}

FECORE_EXPORT  void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 1;
	minor = 0;
	patch = SVNREVISION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  FECoreFactory * PluginGetFactory(int i)
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
	case 17:
		return &plot_matrix_tangent;
	case 18:
		return &plot_matrix_visco_stress;
	case 19:
		return &no_fragment_branching_factory;
	case 20:
		return &psuedo_defered_fragment_branching_factory;
	case 21:
		return &normal_distribution_factory;
	case 22:
		return &default_grow_direction_modifier_factory;
	case 23:
		return &branch_grow_direction_modifier_factory;
	case 24:
		return &gradient_grow_direction_modifier_factory;
	case 25:
		return &anastamosis_grow_direction_modifier_factory;
	case 26:
		return &classic_fragment_seeder_factory;
	case 27:
		return &md_fragment_seeder_factory;
	case 28:
		return &mdbyvol_fragment_seeder_factory;
	case 29:
		return &md_file_seeder_factory;
	case 30:
		return &plot_matrix_elastic_stress;
	case 31:
		return &stopbc_factory;
	case 32:
		return &bouncybc_factory;
	case 33:
		return &same_mbc_factory;
	case 34:
		return &passthrough_mbc_factory;
	case 35:
		return &grow_direction_modifiers_factory;
	case 36:
		return &exponential_distribution_factory;
	case 37:
		return &cauchy_distribution_factory;
	case 38:
		return &chi_squared_distribution_factory;
	default:
		return nullptr;
	}
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginCleanup()
{

}

}
