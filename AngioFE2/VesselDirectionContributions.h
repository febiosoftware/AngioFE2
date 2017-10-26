#pragma once

#include "GrowDirectionModifier.h"

class VesselDirectionContribution : public FEMaterial
{
public:
	explicit VesselDirectionContribution(FEModel * model);
	virtual ~VesselDirectionContribution(){}
	virtual void Update()=0;
	virtual vec3d GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time) =0;
	virtual void SetCulture(Culture * cp);
	double TimescaledFactor(double grow_time) const;
	
protected:
	Culture * culture = nullptr;
	double scale = 1.0;
	bool timedependent;
	
private:
	DECLARE_PARAMETER_LIST();
	
};

class PreviousVesselDirectionContribution : public VesselDirectionContribution
{
public:
	explicit PreviousVesselDirectionContribution(FEModel * model);
	void Update() override{}
	vec3d GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time) override;
};

class FiberVesselDirectionContribution : public VesselDirectionContribution
{
public:
	explicit FiberVesselDirectionContribution(FEModel * model);
	void Update()override{}
	vec3d GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time) override;
};
class ArbitraryVesselDirectionContribution : public VesselDirectionContribution
{
public:
	explicit ArbitraryVesselDirectionContribution(FEModel * model);
	void Update()override;
	void SetCulture(Culture * cp) override;
	vec3d GetContribution(FEAngioMaterialBase* mat, Segment::TIP& tip, double grow_time) override;
private:
	FEPropertyT<GGP> arbitrary_contribution;
};

class VesselDirectionContributionsDirectionModifier : public GrowDirectionModifier
{
public:
	explicit VesselDirectionContributionsDirectionModifier(FEModel * model);
	vec3d GrowModifyGrowDirection(vec3d previous_dir, Segment::TIP& tip, FEAngioMaterialBase* mat, bool branch, double start_time, double grow_time, double& seg_length) override;
	void Update() override;
	void SetCulture(Culture * cp) override;
private:
	FEVecPropertyT<VesselDirectionContribution> vessel_distrubution_contributions;
};