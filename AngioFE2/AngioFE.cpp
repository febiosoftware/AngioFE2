// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

//#define FECORE_DLL
#define FECORE_API

#include "FECore/FECoreKernel.h"
#include "AngioFETask.h"
#include "FECore/FECoreFactory.h"
#include "FESproutBodyForce.h"
#include "FEAngioMaterial.h"
#include "AngioPlot.h"
#include "FEPressureMaterial.h"
#include "VesselDirectionContributions.h"
#ifdef SVN
//#include "svnrev.h"
#define SVNREVISION 0
#else
#define SVNREVISION 0
#endif

//-----------------------------------------------------------------------------
/*
FEPluginFactory_T<AngioFETask       , FETASK_ID    > angiofe_task_factory("angio"   );
FEPluginFactory_T<FESproutBodyForce , FEBODYLOAD_ID> angio_sprout_factory("sprout"  );
FEPluginFactory_T<FEAngioMaterial   , FEMATERIAL_ID> angio_mat_factory   ("angio_mat"   );
FEPluginFactory_T<FEPressureMaterial, FEMATERIAL_ID> pressure_mat_factory("pressure");
FEPluginFactory_T<CommonAngioProperties, FEMATERIAL_ID> common_angio_properties_factory("angio_properties");

FEPluginFactory_T<NoFragmentBranching           , FEMATERIAL_ID> no_fragment_branching_factory            ("no_branch"            );
FEPluginFactory_T<PseudoDeferredFragmentBranching, FEMATERIAL_ID> pseudo_deferred_fragment_branching_factory("pseudo_deferred_branch");

FEPluginFactory_T<FENormalDistribution     , FEMATERIAL_ID> normal_distribution_factory     ("normal_distribution"     );
FEPluginFactory_T<FEUniformDistribution    , FEMATERIAL_ID> uniform_distribution_factory    ("uniform_distribution"    );
FEPluginFactory_T<FEExponentialDistribution, FEMATERIAL_ID> exponential_distribution_factory("exponential_distribution");
FEPluginFactory_T<FECauchyDistribution     , FEMATERIAL_ID> cauchy_distribution_factory     ("cauchy_distribution"     );
FEPluginFactory_T<FEChiSquaredDistribution , FEMATERIAL_ID> chi_squared_distribution_factory("chi_squared_distribution");
FEPluginFactory_T<FEWeibullDistribution    , FEMATERIAL_ID> weibull_distribution_factory    ("weibull_distribution"    );
FEPluginFactory_T<FEGammaDistribution      , FEMATERIAL_ID> gamma_distribution_factory      ("gamma_distribution"      );


FEPluginFactory_T<GrowDirectionModifiers, FEMATERIAL_ID> grow_direction_modifiers_factory("grow_direction_modifiers");

FEPluginFactory_T<DefaultGrowDirectionModifier    , FEMATERIAL_ID> default_grow_direction_modifier_factory    ("default_grow_direction"    );
FEPluginFactory_T<SelectingGrowDirectionModifier, FEMATERIAL_ID> selecting_grow_direction_modifier_factory("selecting_grow_direction");
FEPluginFactory_T<BaseFiberAwareGrowDirectionModifier, FEMATERIAL_ID> base_fiber_grow_direction_modifier_factory("base_fiber_grow_direction");
FEPluginFactory_T<BranchGrowDirectionModifier     , FEMATERIAL_ID> branch_grow_direction_modifier_factory     ("branch_grow_direction"     );
FEPluginFactory_T<RandomBranchGrowDirectionModifier, FEMATERIAL_ID> random_branch_grow_direction_modifier_factory("random_branch_grow_direction");
FEPluginFactory_T<RandomBranchGrowDirectionModifier, FEMATERIAL_ID> random_theta_branch_grow_direction_modifier_factory("random_theta_branch_grow_direction");

FEPluginFactory_T<GradientGrowDirectionModifier   , FEMATERIAL_ID> gradient_grow_direction_modifier_factory   ("density_gradient_grow_direction"   );
FEPluginFactory_T<AnastamosisGrowDirectionModifier, FEMATERIAL_ID> anastamosis_grow_direction_modifier_factory("anastamosis_grow_direction");
FEPluginFactory_T<EdgeDeflectorGrowDirectionModifier, FEMATERIAL_ID> edge_deflector_grow_direction_modifier_factory("edge_deflector_grow_direction");
FEPluginFactory_T<VesselDirectionContributionsDirectionModifier, FEMATERIAL_ID> vessel_direction_contributions_grow_direction_modifier_factory("vessel_direction_contributions");


FEPluginFactory_T<PreviousVesselDirectionContribution, FEMATERIAL_ID> previous_vessel_direction_contribution_factory("previous_vessel_direction");
FEPluginFactory_T<FiberVesselDirectionContribution, FEMATERIAL_ID> fiber_vessel_direction_contribution_factory("fiber_vessel_direction");
FEPluginFactory_T<ArbitraryVesselDirectionContribution, FEMATERIAL_ID> arbitrary_vessel_direction_contribution_factory("arbitrary_vessel_direction");

FEPluginFactory_T<UnitLengthGrowDirectionModifier, FEMATERIAL_ID> unit_length_grow_direction_modifier_factory("unit_length");
FEPluginFactory_T<DensityScaleGrowDirectionModifier, FEMATERIAL_ID> density_scale_grow_direction_modifier_factory("density_length");
FEPluginFactory_T<SegmentLengthGrowDirectionModifier, FEMATERIAL_ID> segment_length_grow_direction_modifier_factory("segment_length");

FEPluginFactory_T<Plot2GGP, FEMATERIAL_ID> plot2_ggp_factory("plot2_ggp");
FEPluginFactory_T<GradientPlot2GGP, FEMATERIAL_ID> gradient_plot2_ggp_factory("gradient_plot2_ggp");
FEPluginFactory_T<MatrixConverterGGP, FEMATERIAL_ID> matrix_converter_ggp_factory("matrix_converter_ggp");
FEPluginFactory_T<ForkedGGP, FEMATERIAL_ID> forked_ggp_factory("forked_ggp");
FEPluginFactory_T<MatrixMixGGP, FEMATERIAL_ID> matrix_mix_ggp_factory("matrix_mix_ggp");
FEPluginFactory_T<EigenValuesGGP, FEMATERIAL_ID> eigen_values_ggp_factory("eigen_values_ggp");
FEPluginFactory_T<EigenVectorsGGP, FEMATERIAL_ID> eigen_vectors_ggp_factory("eigen_vectors_ggp");
FEPluginFactory_T<CrossGGP, FEMATERIAL_ID> cross_ggp_factory("cross_ggp");
FEPluginFactory_T<ThresholdGGP, FEMATERIAL_ID> threshold_ggp_factory("threshold_ggp");
FEPluginFactory_T<ArcCosGGP, FEMATERIAL_ID> arc_cos_ggp_factory("arccos_ggp");
FEPluginFactory_T<ArcSinGGP, FEMATERIAL_ID> arc_sin_ggp_factory("arcsin_ggp");
FEPluginFactory_T<CosGGP, FEMATERIAL_ID> cos_ggp_factory("cos_ggp");
FEPluginFactory_T<SinGGP, FEMATERIAL_ID> sin_ggp_factory("sin_ggp");
FEPluginFactory_T<UnitGGP, FEMATERIAL_ID> unit_diagonal_factory("unit_diagonal_ggp");
FEPluginFactory_T<SetterGGP, FEMATERIAL_ID> setter_ggp_factory("setter_ggp");
FEPluginFactory_T<MatrixSetterGGP, FEMATERIAL_ID> matrix_setter_ggp_factory("matrix_setter_ggp");
FEPluginFactory_T<MatrixInverseGGP, FEMATERIAL_ID> matrix_inverse_ggp_factory("matrix_inverse_ggp");
FEPluginFactory_T<AssertGGP, FEMATERIAL_ID> assert_ggp_factory("assert_ggp");
FEPluginFactory_T<NodalDataGGP, FEMATERIAL_ID> nodal_data_ggp_factory("nodal_data_ggp");
FEPluginFactory_T<NodalDataGradientGGP, FEMATERIAL_ID> nodal_data_gradient_ggp_factory("nodal_data_gradient_ggp");
FEPluginFactory_T<DirectionChangeGGP, FEMATERIAL_ID> direction_change_ggp_factory("direction_change_ggp");


FEPluginFactory_T<ClassicFragmentSeeder      , FEMATERIAL_ID> classic_fragment_seeder_factory("classic"   );
FEPluginFactory_T<MDByVolumeFragmentSeeder   , FEMATERIAL_ID> mdbyvol_fragment_seeder_factory("MDbyVolume");
FEPluginFactory_T<MultiDomainFragmentSeeder  , FEMATERIAL_ID> md_fragment_seeder_factory     ("MD"        );
FEPluginFactory_T<MDAngVessFileFragmentSeeder, FEMATERIAL_ID> md_file_seeder_factory         ("from_file" );

FEPluginFactory_T<BouncyBC, FEMATERIAL_ID> bouncybc_factory("bouncy");
FEPluginFactory_T<StopBC  , FEMATERIAL_ID> stopbc_factory  ("stop"  );

FEPluginFactory_T<SameMBC       , FEMATERIAL_ID> same_mbc_factory       ("same"        );
FEPluginFactory_T<PassThroughMBC, FEMATERIAL_ID> passthrough_mbc_factory("pass_through");

FEPluginFactory_T<FEPlotAngioStress          , FEPLOTDATA_ID> plot_angio_stress            ("angio stress"          );
FEPluginFactory_T<FEPlotVesselStress         , FEPLOTDATA_ID> plot_vessel_stress           ("vessel stress"         );
FEPluginFactory_T<FEPlotMatrixStress         , FEPLOTDATA_ID> plot_matrix_stress           ("matrix stress"         );
FEPluginFactory_T<FEPlotVesselWeight         , FEPLOTDATA_ID> plot_vessel_weight           ("vessel weight"         );
FEPluginFactory_T<FEPlotMatrixWeight         , FEPLOTDATA_ID> plot_matrix_weight           ("matrix weight"         );
FEPluginFactory_T<FEPlotMatrixTangent        , FEPLOTDATA_ID> plot_matrix_tangent          ("matrix tangent"        );
FEPluginFactory_T<FEPlotMatrixViscoStress    , FEPLOTDATA_ID> plot_matrix_visco_stress     ("matrix visco stress"   );
FEPluginFactory_T<FEPlotMatrixElasticStress  , FEPLOTDATA_ID> plot_matrix_elastic_stress   ("matrix elastic stress" );

FEPluginFactory_T<FEPlotAngioECMDensity      , FEPLOTDATA_ID> plot_angio_ecm               ("angio ECM density"     );
FEPluginFactory_T<FEPlotAngioECMAlpha        , FEPLOTDATA_ID> plot_angio_alpha             ("angio ECM alpha"       );
FEPluginFactory_T<FEPlotAngioGradient        , FEPLOTDATA_ID> plot_angio_gradient          ("angio gradient"        );
FEPluginFactory_T<FEPlotAngioGradientCenter  , FEPLOTDATA_ID> plot_angio_gradient_center   ("angio gradient center" );
FEPluginFactory_T<FEPlotAngioMaterialHop     , FEPLOTDATA_ID> plot_angio_material_hop      ("angio material hop"    );
FEPluginFactory_T<FEPlotAngioSegmentBadGrowth, FEPLOTDATA_ID> plot_angio_segment_bad_growth("angio bad growth"      );
FEPluginFactory_T<FEPlotMatrixConectrationGradient, FEPLOTDATA_ID> plot_matrix_concentration_gradient("matrix concecntration gradient");
FEPluginFactory_T<FEPlotMatrixSBMConectration, FEPLOTDATA_ID> plot_matrix_sbm_concentration("matrix sbm concecntration");

FEPluginFactory_T<FEPlotMatrixElastic_m_Q, FEPLOTDATA_ID> plot_matrix_elastic_m_Q("matrix elastic mQ");


FEPluginFactory_T<NullFiberInitializer, FEMATERIAL_ID> null_fiber_initializer("null_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializer, FEMATERIAL_ID> random_fiber_initializer("random_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializerNonMangling, FEMATERIAL_ID> random_fiber_initializer_non_mangling("random_fiber_initializer_non_mangling");
FEPluginFactory_T<ExplicitDistributionsFiberInitializer, FEMATERIAL_ID> explicit_distribution_fiber_initializer("explicit_distribution_fiber_initializer");
FEPluginFactory_T<RandomFiberInitializerPE, FEMATERIAL_ID> random_fiber_initializer_pe("random_fiber_initializer_pe");
FEPluginFactory_T<EllipsoidalFiberInitializer, FEMATERIAL_ID> ellipsoidal_fiber_initializer("ellipsoidal_fiber_initializer");
*/
//-----------------------------------------------------------------------------
FECORE_EXPORT  unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	REGISTER_FECORE_CLASS(AngioFETask, FETASK_ID, "angio");
	REGISTER_FECORE_CLASS(FESproutBodyForce, FEBODYLOAD_ID, "sprout");
	REGISTER_FECORE_CLASS(FEAngioMaterial, FEMATERIAL_ID, "angio_mat");
	REGISTER_FECORE_CLASS(FEPressureMaterial, FEMATERIAL_ID, "pressure");
	REGISTER_FECORE_CLASS(CommonAngioProperties, FEMATERIAL_ID, "angio_properties");

	REGISTER_FECORE_CLASS(NoFragmentBranching, FEMATERIAL_ID, "no_branch");
	REGISTER_FECORE_CLASS(PseudoDeferredFragmentBranching, FEMATERIAL_ID, "pseudo_deferred_branch");

	REGISTER_FECORE_CLASS(FENormalDistribution, FEMATERIAL_ID, "normal_distribution");
	REGISTER_FECORE_CLASS(FEUniformDistribution, FEMATERIAL_ID, "uniform_distribution");
	REGISTER_FECORE_CLASS(FEExponentialDistribution, FEMATERIAL_ID, "exponential_distribution");
	REGISTER_FECORE_CLASS(FECauchyDistribution, FEMATERIAL_ID, "cauchy_distribution");
	REGISTER_FECORE_CLASS(FEChiSquaredDistribution, FEMATERIAL_ID, "chi_squared_distribution");
	REGISTER_FECORE_CLASS(FEWeibullDistribution, FEMATERIAL_ID, "weibull_distribution");
	REGISTER_FECORE_CLASS(FEGammaDistribution, FEMATERIAL_ID, "gamma_distribution");


	REGISTER_FECORE_CLASS(GrowDirectionModifiers, FEMATERIAL_ID, "grow_direction_modifiers");

	REGISTER_FECORE_CLASS(DefaultGrowDirectionModifier, FEMATERIAL_ID, "default_grow_direction");
	REGISTER_FECORE_CLASS(SelectingGrowDirectionModifier, FEMATERIAL_ID, "selecting_grow_direction");
	REGISTER_FECORE_CLASS(BaseFiberAwareGrowDirectionModifier, FEMATERIAL_ID, "base_fiber_grow_direction");
	REGISTER_FECORE_CLASS(BranchGrowDirectionModifier, FEMATERIAL_ID, "branch_grow_direction");
	REGISTER_FECORE_CLASS(RandomBranchGrowDirectionModifier, FEMATERIAL_ID, "random_branch_grow_direction");
	REGISTER_FECORE_CLASS(RandomThetaBranchGrowDirectionModifier, FEMATERIAL_ID, "random_theta_branch_grow_direction");

	REGISTER_FECORE_CLASS(GradientGrowDirectionModifier, FEMATERIAL_ID, "density_gradient_grow_direction");
	REGISTER_FECORE_CLASS(AnastamosisGrowDirectionModifier, FEMATERIAL_ID, "anastamosis_grow_direction");
	REGISTER_FECORE_CLASS(EdgeDeflectorGrowDirectionModifier, FEMATERIAL_ID, "edge_deflector_grow_direction");
	REGISTER_FECORE_CLASS(VesselDirectionContributionsDirectionModifier, FEMATERIAL_ID, "vessel_direction_contributions");


	REGISTER_FECORE_CLASS(PreviousVesselDirectionContribution, FEMATERIAL_ID, "previous_vessel_direction");
	REGISTER_FECORE_CLASS(FiberVesselDirectionContribution, FEMATERIAL_ID, "fiber_vessel_direction");
	REGISTER_FECORE_CLASS(ArbitraryVesselDirectionContribution, FEMATERIAL_ID, "arbitrary_vessel_direction");

	REGISTER_FECORE_CLASS(UnitLengthGrowDirectionModifier, FEMATERIAL_ID, "unit_length");
	REGISTER_FECORE_CLASS(DensityScaleGrowDirectionModifier, FEMATERIAL_ID, "density_length");
	REGISTER_FECORE_CLASS(RefDensityScaleGrowDirectionModifier, FEMATERIAL_ID, "ref_density_length");
	REGISTER_FECORE_CLASS(SegmentLengthGrowDirectionModifier, FEMATERIAL_ID, "segment_length");

	REGISTER_FECORE_CLASS(Plot2GGP, FEMATERIAL_ID, "plot2_ggp");
	REGISTER_FECORE_CLASS(GradientPlot2GGP, FEMATERIAL_ID, "gradient_plot2_ggp");
	REGISTER_FECORE_CLASS(MatrixConverterGGP, FEMATERIAL_ID, "matrix_converter_ggp");
	REGISTER_FECORE_CLASS(ForkedGGP, FEMATERIAL_ID, "forked_ggp");
	REGISTER_FECORE_CLASS(MatrixMixGGP, FEMATERIAL_ID, "matrix_mix_ggp");
	REGISTER_FECORE_CLASS(EigenValuesGGP, FEMATERIAL_ID, "eigen_values_ggp");
	REGISTER_FECORE_CLASS(EigenVectorsGGP, FEMATERIAL_ID, "eigen_vectors_ggp");
	REGISTER_FECORE_CLASS(CrossGGP, FEMATERIAL_ID, "cross_ggp");
	REGISTER_FECORE_CLASS(ThresholdGGP, FEMATERIAL_ID, "threshold_ggp");
	REGISTER_FECORE_CLASS(ArcCosGGP, FEMATERIAL_ID, "arccos_ggp");
	REGISTER_FECORE_CLASS(ArcSinGGP, FEMATERIAL_ID, "arcsin_ggp");
	REGISTER_FECORE_CLASS(CosGGP, FEMATERIAL_ID, "cos_ggp");
	REGISTER_FECORE_CLASS(SinGGP, FEMATERIAL_ID, "sin_ggp");
	REGISTER_FECORE_CLASS(UnitGGP, FEMATERIAL_ID, "unit_diagonal_ggp");
	REGISTER_FECORE_CLASS(SetterGGP, FEMATERIAL_ID, "setter_ggp");
	REGISTER_FECORE_CLASS(MatrixSetterGGP, FEMATERIAL_ID, "matrix_setter_ggp");
	REGISTER_FECORE_CLASS(MatrixInverseGGP, FEMATERIAL_ID, "matrix_inverse_ggp");
	REGISTER_FECORE_CLASS(AssertGGP, FEMATERIAL_ID, "assert_ggp");
	REGISTER_FECORE_CLASS(NodalDataGGP, FEMATERIAL_ID, "nodal_data_ggp");
	REGISTER_FECORE_CLASS(NodalDataGradientGGP, FEMATERIAL_ID, "nodal_data_gradient_ggp");
	REGISTER_FECORE_CLASS(DirectionChangeGGP, FEMATERIAL_ID, "direction_change_ggp");


	REGISTER_FECORE_CLASS(ClassicFragmentSeeder, FEMATERIAL_ID, "classic");
	REGISTER_FECORE_CLASS(MDByVolumeFragmentSeeder, FEMATERIAL_ID, "MDbyVolume");
	REGISTER_FECORE_CLASS(MultiDomainFragmentSeeder, FEMATERIAL_ID, "MD");
	REGISTER_FECORE_CLASS(MDAngVessFileFragmentSeeder, FEMATERIAL_ID, "from_file");

	REGISTER_FECORE_CLASS(BouncyBC, FEMATERIAL_ID, "bouncy");
	REGISTER_FECORE_CLASS(StopBC, FEMATERIAL_ID, "stop");

	REGISTER_FECORE_CLASS(SameMBC, FEMATERIAL_ID, "same");
	REGISTER_FECORE_CLASS(PassThroughMBC, FEMATERIAL_ID, "pass_through");

	REGISTER_FECORE_CLASS(FEPlotAngioStress, FEPLOTDATA_ID, "angio stress");
	REGISTER_FECORE_CLASS(FEPlotVesselStress, FEPLOTDATA_ID, "vessel stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixStress, FEPLOTDATA_ID, "matrix stress");
	REGISTER_FECORE_CLASS(FEPlotVesselWeight, FEPLOTDATA_ID, "vessel weight");
	REGISTER_FECORE_CLASS(FEPlotMatrixWeight, FEPLOTDATA_ID, "matrix weight");
	REGISTER_FECORE_CLASS(FEPlotMatrixTangent, FEPLOTDATA_ID,"matrix tangent");
	REGISTER_FECORE_CLASS(FEPlotMatrixViscoStress, FEPLOTDATA_ID, "matrix visco stress");
	REGISTER_FECORE_CLASS(FEPlotMatrixElasticStress, FEPLOTDATA_ID, "matrix elastic stress");

	REGISTER_FECORE_CLASS(FEPlotAngioECMDensity, FEPLOTDATA_ID, "angio ECM density");
	REGISTER_FECORE_CLASS(FEPlotAngioECMAlpha, FEPLOTDATA_ID, "angio ECM alpha");
	REGISTER_FECORE_CLASS(FEPlotAngioGradient, FEPLOTDATA_ID, "angio gradient");
	REGISTER_FECORE_CLASS(FEPlotAngioGradientCenter, FEPLOTDATA_ID, "angio gradient center");
	REGISTER_FECORE_CLASS(FEPlotAngioMaterialHop, FEPLOTDATA_ID, "angio material hop");
	REGISTER_FECORE_CLASS(FEPlotAngioSegmentBadGrowth, FEPLOTDATA_ID, "angio bad growth");
	REGISTER_FECORE_CLASS(FEPlotMatrixConectrationGradient, FEPLOTDATA_ID, "matrix concecntration gradient");
	REGISTER_FECORE_CLASS(FEPlotMatrixSBMConectration, FEPLOTDATA_ID, "matrix sbm concecntration");

	REGISTER_FECORE_CLASS(FEPlotMatrixElastic_m_Q, FEPLOTDATA_ID, "matrix elastic mQ");


	REGISTER_FECORE_CLASS(NullFiberInitializer, FEMATERIAL_ID, "null_fiber_initializer");
	REGISTER_FECORE_CLASS(RandomFiberInitializer, FEMATERIAL_ID, "random_fiber_initializer");
	REGISTER_FECORE_CLASS(RandomFiberInitializerNonMangling, FEMATERIAL_ID, "random_fiber_initializer_non_mangling");
	REGISTER_FECORE_CLASS(ExplicitDistributionsFiberInitializer, FEMATERIAL_ID, "explicit_distribution_fiber_initializer");
	REGISTER_FECORE_CLASS(RandomFiberInitializerPE, FEMATERIAL_ID, "random_fiber_initializer_pe");
	REGISTER_FECORE_CLASS(EllipsoidalFiberInitializer, FEMATERIAL_ID, "ellipsoidal_fiber_initializer");
}

FECORE_EXPORT  void GetPluginVersion(int & major, int & minor, int & patch)
{
	major = 2;
	minor = 1;
	patch = SVNREVISION;
}

//-----------------------------------------------------------------------------
/*FECORE_EXPORT  FECoreFactory * PluginGetFactory(int i)
{
	std::vector<FECoreFactory *> addon_classes{ 
		&angiofe_task_factory, &angio_sprout_factory, &angio_mat_factory,
		&angio_mat_factory, &pressure_mat_factory,
		//plot classes
		&plot_angio_stress, &plot_angio_stress,
		&plot_angio_ecm, &plot_angio_alpha, &plot_angio_gradient, &plot_angio_gradient_center,
		&plot_angio_material_hop, &plot_angio_segment_bad_growth, &plot_vessel_stress, &plot_matrix_stress,
		&plot_vessel_weight, &plot_matrix_weight, &plot_matrix_tangent, &plot_matrix_visco_stress,
		&plot_matrix_elastic_m_Q, &plot_matrix_elastic_stress,

		//stopgap plot classes
		&plot_matrix_concentration_gradient, &plot_matrix_sbm_concentration,

		//fiber initializers
		&null_fiber_initializer, &random_fiber_initializer,
		&random_fiber_initializer_non_mangling, &explicit_distribution_fiber_initializer,
		&random_fiber_initializer_pe, &ellipsoidal_fiber_initializer,
		
		//branching factories
		&no_fragment_branching_factory, &pseudo_deferred_fragment_branching_factory,
		//grow direction modifiers
		&grow_direction_modifiers_factory, &base_fiber_grow_direction_modifier_factory,
		&unit_length_grow_direction_modifier_factory, &segment_length_grow_direction_modifier_factory,
		&default_grow_direction_modifier_factory, &selecting_grow_direction_modifier_factory,
		&branch_grow_direction_modifier_factory,
		&random_branch_grow_direction_modifier_factory, &random_theta_branch_grow_direction_modifier_factory,
		&gradient_grow_direction_modifier_factory, &anastamosis_grow_direction_modifier_factory,
		&density_scale_grow_direction_modifier_factory, &edge_deflector_grow_direction_modifier_factory,
		&vessel_direction_contributions_grow_direction_modifier_factory,

		//vessel contribution modifiers
		&previous_vessel_direction_contribution_factory, &fiber_vessel_direction_contribution_factory,
		&arbitrary_vessel_direction_contribution_factory,

		//fragment seeders
		&classic_fragment_seeder_factory,  &md_fragment_seeder_factory, &mdbyvol_fragment_seeder_factory,
		&md_file_seeder_factory,
		//boundary conditions
		&stopbc_factory, &bouncybc_factory, &same_mbc_factory,
		&passthrough_mbc_factory, 
		//random distribution
		&cauchy_distribution_factory, &chi_squared_distribution_factory, &weibull_distribution_factory,
		&gamma_distribution_factory,&normal_distribution_factory, &exponential_distribution_factory,
		&uniform_distribution_factory,
		//ggp's
		&plot2_ggp_factory, &gradient_plot2_ggp_factory,
		&matrix_converter_ggp_factory, &forked_ggp_factory, &cross_ggp_factory,
		&threshold_ggp_factory, &nodal_data_ggp_factory, &nodal_data_gradient_ggp_factory,
		&arc_cos_ggp_factory, &arc_sin_ggp_factory, &cos_ggp_factory, &sin_ggp_factory,
		&matrix_inverse_ggp_factory, &eigen_vectors_ggp_factory, &eigen_values_ggp_factory,
		&setter_ggp_factory, &matrix_setter_ggp_factory,
		&assert_ggp_factory, &unit_diagonal_factory, &direction_change_ggp_factory,
		&matrix_mix_ggp_factory,

		//other needed items
		&common_angio_properties_factory
	};

	if(i < addon_classes.size())
	{
		return addon_classes[i];
	}
	return nullptr;

}
*/
//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginCleanup()
{

}

