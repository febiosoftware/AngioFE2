#pragma once
#include "StdAfx.h"
#include <FECore/FEDomain.h>
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
class FragmentSeeder
{
public:
	virtual bool SeedFragments(SimulationTime& time, Culture * culture) = 0;
	FragmentSeeder(CultureParameters * cp, FEAngio & angio) : m_angio(angio), culture_params(cp) {}
	virtual ~FragmentSeeder(){}
protected:
	FEAngio & m_angio;
	CultureParameters * culture_params;
	std::vector<FEDomain *> domains;
	virtual bool createInitFrag(Segment& Seg, SegGenItem & item, Culture * culture);
};

class ClassicFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	ClassicFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
	// Seed an initial fragment within the grid
	bool createInitFrag(Segment& Seg);
};
class MultiDomainFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MultiDomainFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
};

class MDByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDByVolumeFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
};

class MDAngVessFileFragmentSeeder :public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDAngVessFileFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
	std::ifstream infile;
};