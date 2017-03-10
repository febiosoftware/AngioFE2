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
	FEParam * m = GetParameter(ParamString("mean"));
	FEParam * s = GetParameter(ParamString("stddev"));

	int mlci = m->GetLoadCurve();
	int slci = s->GetLoadCurve();

	FEModel * model = GetFEModel();
	if (mlci > 0)
	{
		FEDataLoadCurve * mlc = dynamic_cast<FEDataLoadCurve*>(model->GetLoadCurve(mlci));
		mlc->SetInterpolation(FEDataLoadCurve::STEP);
	}
	if (slci > 0)
	{
		FEDataLoadCurve * slc = dynamic_cast<FEDataLoadCurve*>( model->GetLoadCurve(slci));
		slc->SetInterpolation(FEDataLoadCurve::STEP);
	}
	return true;
}

void FENormalDistribution::StepToTime(double time)
{
	FEParam * m = GetParameter(ParamString("mean"));
	FEParam * s = GetParameter(ParamString("stddev"));

	int mlci = m->GetLoadCurve();
	int slci = s->GetLoadCurve();

	bool change = false;
	FEModel * model = GetFEModel();
	double mv = mean;
	double sv = stddev;
	if (mlci > 0)
	{
		FELoadCurve * mlc = model->GetLoadCurve(mlci);
		mv =  mlc->Value(time);
		if (mv != prev_mean)
		{
			change = true;
		}
	}
	if (slci > 0)
	{
		FELoadCurve * slc = model->GetLoadCurve(slci);
		sv = slc->Value(time);
		if (sv != prev_stddev)
		{
			change = true;
		}
	}
	if (change)
	{
		//rebuild the distribution
		prev_mean = mv;
		prev_stddev = sv;
		nd = std::normal_distribution<double>(mv, sv);
	}
}

BEGIN_PARAMETER_LIST(FEProbabilityDistribution, FEMaterial)
ADD_PARAMETER(max_retries, FE_PARAM_INT, "max_retries");
END_PARAMETER_LIST();

BEGIN_PARAMETER_LIST(FENormalDistribution, FEProbabilityDistribution)
ADD_PARAMETER(mean, FE_PARAM_DOUBLE, "mean");
ADD_PARAMETER(stddev, FE_PARAM_DOUBLE, "stddev");

END_PARAMETER_LIST();