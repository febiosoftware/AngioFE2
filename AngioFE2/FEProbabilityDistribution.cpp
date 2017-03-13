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
	FEParam * m = GetParameter(ParamString("mean"));
	FEParam * s = GetParameter(ParamString("stddev"));

	int mlci = m->GetLoadCurve();
	int slci = s->GetLoadCurve();

	bool change = ChangeInParam("mean", time, prev_mean, mean) || ChangeInParam("stddev", time, prev_stddev, stddev);
	if (change)
	{
		//rebuild the distribution
		prev_mean = mean;
		prev_stddev = stddev;
		nd = std::normal_distribution<double>(mean, stddev);
	}
}

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

BEGIN_PARAMETER_LIST(FENormalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(mean, FE_PARAM_DOUBLE, "mean");
ADD_PARAMETER(stddev, FE_PARAM_DOUBLE, "stddev");
END_PARAMETER_LIST();



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
	return true;
}

void FEExponentialDistribution::StepToTime(double time)
{
	bool change = ChangeInParam("lambda", time, prev_lambda, lambda);
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