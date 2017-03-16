#include "FEProbabilityDistribution.h"
#include "FECore/FEModel.h"
#include "FECore/LoadCurve.h"
#include <FECore/FEDataLoadCurve.h>

//implemenations of FENormalDistribution
double FENormalDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = nd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FENormalDistribution::Init()
{
	nd = std::normal_distribution<double>(mean, stddev);
	prev_mean = mean;
	prev_stddev = stddev;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("mean");
	SetLoadCurveToStep("stddev");

	return true;
}

void FENormalDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("mean", time, prev_mean, mean) || ChangeInParam("stddev", time, prev_stddev, stddev);
	if (change)
	{
		//rebuild the distribution
		prev_mean = mean;
		prev_stddev = stddev;
		nd = std::normal_distribution<double>(mean, stddev);
	}
}


BEGIN_PARAMETER_LIST(FENormalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(mean, FE_PARAM_DOUBLE, "mean");
ADD_PARAMETER(stddev, FE_PARAM_DOUBLE, "stddev");
END_PARAMETER_LIST();

BEGIN_PARAMETER_LIST(FEProbabilityDistribution, FEMaterial)
ADD_PARAMETER(max_retries, FE_PARAM_INT, "max_retries");
END_PARAMETER_LIST();

void FEProbabilityDistribution::SetLoadCurveToStep(const char * param)
{
	//if load curves are used they must use step interpolation
	FEParam * m = GetParameter(ParamString(param));

	int mlci = m->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci > 0)
	{
		FEDataLoadCurve * mlc = dynamic_cast<FEDataLoadCurve*>(model->GetLoadCurve(mlci));
		mlc->SetInterpolation(FEDataLoadCurve::STEP);
	}
}

bool FEProbabilityDistribution::ChangeInParam(const char * param, double time, double & prev, double & new_p)
{
	FEParam * m = GetParameter(ParamString(param));

	int mlci = m->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci > 0)
	{
		FELoadCurve * mlc = model->GetLoadCurve(mlci);
		new_p = mlc->Value(time);
		if (new_p != prev)
		{
			return true;
		}
	}
	return false;
}




//implemenations of FENormalDistribution
double FEExponentialDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = ed(re);
		if (val > 0.0)
			return mult*val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FEExponentialDistribution::Init()
{
	ed = std::exponential_distribution<double>(lambda);
	prev_lambda = lambda;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("lambda");
	SetLoadCurveToStep("mult");
	return true;
}

void FEExponentialDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("lambda", time, prev_lambda, lambda);
	ChangeInParam("mult", time, prev_mult, mult);
	if (change)
	{
		//rebuild the distribution
		prev_lambda = lambda;
		ed = std::exponential_distribution<double>(lambda);
	}
}

BEGIN_PARAMETER_LIST(FEExponentialDistribution, FEProbabilityDistribution)
ADD_PARAMETER(lambda, FE_PARAM_DOUBLE, "lambda");
ADD_PARAMETER(mult, FE_PARAM_DOUBLE, "mult");
END_PARAMETER_LIST();

//implemenations of FENormalDistribution
double FECauchyDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = cd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FECauchyDistribution::Init()
{
	cd = std::cauchy_distribution<double>(a, b);
	prev_a = a;
	prev_b = b;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");

	return true;
}

void FECauchyDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("a", time, prev_a, a) || ChangeInParam("b", time, prev_b, b);
	if (change)
	{
		//rebuild the distribution
		prev_a = a;
		prev_b = b;
		cd = std::cauchy_distribution<double>(a, b);
	}
}


BEGIN_PARAMETER_LIST(FECauchyDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();


//implemenations of FENormalDistribution
double FEChiSquaredDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = cd(re);
		if (val > 0.0)
			return mult*val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FEChiSquaredDistribution::Init()
{
	cd = std::chi_squared_distribution<double>(dof);
	prev_dof = dof;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("dof");
	SetLoadCurveToStep("mult");

	return true;
}

void FEChiSquaredDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("dof", time, prev_dof, dof);
	ChangeInParam("mult", time, prev_mult, mult);
	if (change)
	{
		//rebuild the distribution
		prev_dof = dof;
		cd = std::chi_squared_distribution<double>(dof);
	}
}


BEGIN_PARAMETER_LIST(FEChiSquaredDistribution, FEProbabilityDistribution)
ADD_PARAMETER(dof, FE_PARAM_DOUBLE, "dof");
ADD_PARAMETER(mult, FE_PARAM_DOUBLE, "mult");
END_PARAMETER_LIST();


//implemenations of FENormalDistribution
double FEWeibullDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = wd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FEWeibullDistribution::Init()
{
	wd = std::weibull_distribution<double>(a, b);
	prev_a = a;
	prev_b = b;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("a");
	SetLoadCurveToStep("b");

	return true;
}

void FEWeibullDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("a", time, prev_a, a) || ChangeInParam("b", time, prev_b, b);
	if (change)
	{
		//rebuild the distribution
		prev_a = a;
		prev_b = b;
		wd = std::weibull_distribution<double>(a, b);
	}
}


BEGIN_PARAMETER_LIST(FEWeibullDistribution, FEProbabilityDistribution)
ADD_PARAMETER(a, FE_PARAM_DOUBLE, "a");
ADD_PARAMETER(b, FE_PARAM_DOUBLE, "b");
END_PARAMETER_LIST();


//implemenations of FENormalDistribution
double FEGammaDistribution::NextValue(angiofe_random_engine & re)
{
	for (int i = 0; i < max_retries; i++)
	{
		double val = gd(re);
		if (val > 0.0)
			return val;
	}
	return std::numeric_limits<double>::quiet_NaN();
}

bool FEGammaDistribution::Init()
{
	gd = std::gamma_distribution<double>(alpha, beta);
	prev_alpha = alpha;
	prev_beta = beta;
	//if load curves are used they must use step interpolation
	SetLoadCurveToStep("alpha");
	SetLoadCurveToStep("beta");

	return true;
}

void FEGammaDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("alpha", time, prev_alpha, alpha) || ChangeInParam("beta", time, prev_beta, beta);
	if (change)
	{
		//rebuild the distribution
		prev_alpha = alpha;
		prev_beta = beta;
		gd = std::gamma_distribution<double>(alpha, beta);
	}
}


BEGIN_PARAMETER_LIST(FEGammaDistribution, FEProbabilityDistribution)
ADD_PARAMETER(alpha, FE_PARAM_DOUBLE, "alpha");
ADD_PARAMETER(beta, FE_PARAM_DOUBLE, "beta");
END_PARAMETER_LIST();