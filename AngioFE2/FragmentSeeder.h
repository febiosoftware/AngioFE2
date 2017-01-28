#pragma once
#include "StdAfx.h"
#include <FECore/FEDomain.h>
#include <FECore/FEMaterial.h>
#include "CultureParameters.h"


//this header contains the definitions of all of the fragment seeder objects
class SimulationTime;
class Culture;
class FEAngio;
class Segment;

//used for generating initial segments
struct SegGenItem
{
	bool operator<(SegGenItem const & item)const{ return weight < item.weight; }
	double weight;
	FEDomain * domain;
	int ielement;
};
class Culture;
class FragmentSeeder : public FEMaterial
{
public:
	virtual bool SeedFragments(SimulationTime& time, Culture * culture) = 0;
	FragmentSeeder(FEModel * model);
	virtual ~FragmentSeeder(){}
	//must be called before anything else is done but construction
	void SetCulture(Culture * cp);
protected:
	Culture * culture;
	std::vector<FEDomain *> domains;
	virtual bool createInitFrag(Segment& Seg, SegGenItem & item, Culture * culture);
	int number_fragments =0;
	double initial_vessel_length = 20.0;
	DECLARE_PARAMETER_LIST();
};

class ClassicFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	ClassicFragmentSeeder(FEModel * model);
private:
	// Seed an initial fragment within the grid
	bool createInitFrag(Segment& Seg);
};
class MultiDomainFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MultiDomainFragmentSeeder(FEModel * model);
private:
};

class MDByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDByVolumeFragmentSeeder(FEModel * model);
private:
};

class MDAngVessFileFragmentSeeder :public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDAngVessFileFragmentSeeder(FEModel * model);
private:
	DECLARE_PARAMETER_LIST();
	std::ifstream infile;
	char vessel_file[256];
};