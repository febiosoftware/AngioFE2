#pragma once
#include "StdAfx.h"
#include <FECore/FEMaterial.h>
//pure virtual base class for probability distributions
class FEProbabilityDistribution : public FEMaterial
{
public:
	explicit FEProbabilityDistribution(FEModel* pfem) : FEMaterial(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	virtual double NextValue(angiofe_random_engine & re) = 0;

	
	//when the time changes if the distribution needs modified this is the time to do it
	virtual void StepToTime(double time) = 0;

	
protected:
	int max_retries = 10;
	DECLARE_PARAMETER_LIST();
	void SetLoadCurveToStep(const char * param);
	bool ChangeInParam(const char * param, double time, double & prev, double & new_p);
private:
	
};

class FENormalDistribution : public FEProbabilityDistribution
{
public:
	explicit FENormalDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double mean = 1.0;//distribution's mean
	double stddev = 1.0;//distribution's standard deviation

	std::normal_distribution<double> nd;

	double prev_mean = 1.0;
	double prev_stddev = 1.0;
};

class FEUniformDistribution : public FEProbabilityDistribution
{
public:
	explicit FEUniformDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double a = 0.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation

	std::uniform_real_distribution<double> rd;
	bool time_clamped = true;

	double prev_a = 0.0;
	double prev_b = 1.0;
};

class FEExponentialDistribution : public FEProbabilityDistribution
{
public:
	explicit FEExponentialDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double lambda = 1.0;//distribution's lambda
	double mult = 1.0;//multiplier

	std::exponential_distribution<double> ed;

	double prev_lambda = 1.0;
	double prev_mult = mult;
};

class FECauchyDistribution : public FEProbabilityDistribution
{
public:
	explicit FECauchyDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double a = 1.0;//distribution's mean
	double b = 1.0;//distribution's standard deviation

	std::cauchy_distribution<double> cd;

	double prev_a = a;
	double prev_b = b;
};


class FEChiSquaredDistribution : public FEProbabilityDistribution
{
public:
	explicit FEChiSquaredDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double dof = 1.0;//distribution's x^2
	double mult = 1.0;

	std::chi_squared_distribution<double> cd;

	double prev_dof = dof;
	double prev_mult = mult;
};

class FEWeibullDistribution : public FEProbabilityDistribution
{
public:
	explicit FEWeibullDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;

protected:
	DECLARE_PARAMETER_LIST();
private:
	double a = 1.0;
	double b = 1.0;

	std::weibull_distribution<double> wd;

	double prev_a = a;
	double prev_b = b;

	
};

class FEGammaDistribution : public FEProbabilityDistribution
{
public:
	explicit FEGammaDistribution(FEModel* pfem) : FEProbabilityDistribution(pfem) {}

	//generates the next value in the given sequence which fits a given distribution
	//this value cannot be zero or less if the value is zero or less the result will be redrawn up to max_retries
	//nan will be returned if the distribution fails to find a suitable number
	double NextValue(angiofe_random_engine & re) override;

	bool Init() override;

	void StepToTime(double time) override;
protected:
	DECLARE_PARAMETER_LIST();
private:
	double alpha = 1.0;
	double beta = 1.0;

	std::gamma_distribution<double> gd;

	double prev_alpha = alpha;
	double prev_beta = beta;
};