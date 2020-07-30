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
FECORE_EXPORT  void PluginCleanup()
{

}

