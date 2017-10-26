#include "VesselDirectionContributions.h"
#include "Culture.h"
#include "FEAngioMaterialBase.h"


VesselDirectionContributionsDirectionModifier::VesselDirectionContributionsDirectionModifier(FEModel * model): GrowDirectionModifier(model)
{
	AddProperty(&vessel_distrubution_contributions, "vdm");
}

vec3d VesselDirectionContributionsDirectionModifier::GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterialBase* mat, bool branch, double start_time, double grow_time, double& seg_length)
{
	vec3d sum(0, 0, 0);
	for(int i =0; i < vessel_distrubution_contributions.size();i++)
	{
		sum += vessel_distrubution_contributions[i]->GetContribution(mat,tip, grow_time);
	}
	sum.unit();
	return sum;
}
void VesselDirectionContributionsDirectionModifier::SetCulture(Culture * cp)
{
	for (int i = 0; i < vessel_distrubution_contributions.size(); i++)
	{
		vessel_distrubution_contributions[i]->SetCulture(cp);
	}
}

void VesselDirectionContributionsDirectionModifier::Update()
{
	for (int i = 0; i < vessel_distrubution_contributions.size(); i++)
	{
		vessel_distrubution_contributions[i]->Update();
	}
}

VesselDirectionContribution::VesselDirectionContribution(FEModel * model) : FEMaterial(model)
{
	
}

void VesselDirectionContribution::SetCulture(Culture * cp)
{
	culture = cp;
}

vec3d FiberVesselDirectionContribution::GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time)
{
	double lambda;
	return  culture->m_pmat->fiber_manager->GetFiberDirection(tip.pt, lambda) * TimescaledFactor(grow_time);
}

double VesselDirectionContribution::TimescaledFactor(double grow_time) const
{
	if(timedependent)
	{
		return scale*grow_time;
	}
	return scale;
}

FiberVesselDirectionContribution::FiberVesselDirectionContribution(FEModel * model): VesselDirectionContribution(model)
{
	timedependent = true;
}

PreviousVesselDirectionContribution::PreviousVesselDirectionContribution(FEModel * model) : VesselDirectionContribution(model)
{
	timedependent = false;
}
vec3d PreviousVesselDirectionContribution::GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time)
{
	return tip.u * TimescaledFactor(grow_time);
}

BEGIN_PARAMETER_LIST(VesselDirectionContribution, FEMaterial)
ADD_PARAMETER2(scale, FE_PARAM_DOUBLE, FE_RANGE_GREATER_OR_EQUAL(0), "scale");
ADD_PARAMETER(timedependent, FE_PARAM_BOOL, "timedependent");
END_PARAMETER_LIST();


ArbitraryVesselDirectionContribution::ArbitraryVesselDirectionContribution(FEModel * model) : VesselDirectionContribution(model)
{
	timedependent = false;
	AddProperty(&arbitrary_contribution, "arbitrary_contributions");
	arbitrary_contribution.m_brequired = false;
}

void ArbitraryVesselDirectionContribution::Update()
{
	if(arbitrary_contribution)
	{
		arbitrary_contribution->Update();
	}
}

void ArbitraryVesselDirectionContribution::SetCulture(Culture * cp)
{
	if (arbitrary_contribution)
		arbitrary_contribution->SetCulture(cp);
	VesselDirectionContribution::SetCulture(cp);
}

vec3d ArbitraryVesselDirectionContribution::GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time)
{
	if(arbitrary_contribution)
	{
		mat3d id(1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0);
		vec3d in(1, 1, 1);
		mat3d res = arbitrary_contribution->Operation(id, in, mat, tip);
		return res*in*TimescaledFactor(grow_time);
	}
	return vec3d(0,0,0);
}