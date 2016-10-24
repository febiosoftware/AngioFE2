#pragma once

#include <list>
#include "BC.h"
#include "Segment.h"
#include <vector>
using namespace std;

//-----------------------------------------------------------------------------
class FEAngio;
class SimulationTime;
class GridPoint;
class FEAngioMaterial;
class BouncyBC;

//-----------------------------------------------------------------------------
typedef list<Segment> SegmentList;
typedef list<Segment>::iterator SegIter;
typedef list<Segment>::const_iterator ConstSegIter;
typedef list<Segment::TIP*>	SegmentTipList;
typedef list<Segment::TIP*>::iterator TipIter;
typedef list<Segment::TIP*>::const_iterator ConstTipIter;

const int MAXPARAMSIZE = 256;

struct CultureParameters
{
	CultureParameters();                                                         

	//not currently exposed
	//These parameters govern the amount of growth at each time step for the microvessels
	//          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double	m_y0 = -19.1278;   // Bottom of sigmoid curve
	double	m_culture_a = 1900.0;    // Distance between top and bottom of the curve (a + y0 = top of curve)
	double	m_x0 = 4.9474;   // Time point at which t is halfway between top & bottom
	double	m_culture_b = 1.4549;    // Steepness of the curve

	//the initial material parameters a,b,N
	double sprout_s_mag = 3.72e-12;
	double sprout_s_range = 0.004;
	double sprout_s_width = 2.0;

	//dont expose
	double	m_initial_vessel_length = 0.0;	// initial vessel length
	double	m_initial_branch_probability = 0.2697;     // Probability that initial segments will branch
	bool m_bzfibflat = false;//flatten the fibers in the z direction

	//parameters adjustable by the user for vessel network morphology
	// Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjustment = 1.0;
	// If the shortest distance vector between a segment and any other segment is less than this value
	//TODO: I think making this a percentage of the growth length makes more sense.
	double m_anastomosis_distance = 25;
	double	m_branch_chance = 0.1;    // Probability of forming a new branch
	bool	m_branch= false; //whether the vessels are allowed to branch
	bool	m_anastomosis = false; //whether the vessels are allowed to fuse together
	// vessel_width - Diameter of microvessels (Default: 7 um)
	double m_vessel_width = 7.0;
	char m_boundary_condition_type[MAXPARAMSIZE];//currently s for stop, or b for bouncy 

	//parameters for ECM density/alignment
	int m_matrix_condition = 0; // flag indicating how the collagen fibers are oriented initially ( 0 = random classic mode, 1 multimaterial mode,  3 = along local element direction)
	double m_matrix_density = 3.0;//mg/ml 3.0 default

	//symmety
	vec3d m_symmetry_plane = vec3d(0.0, 0.0, 0.0); //symmetry plane

	//uncategoriezed
	int m_composite_material = 1;//whether or not the material is a composite material
	
	//TODO: doesn't appear to be the same as sprout_s_mag
	double m_sprout_force = 1.0;

	int m_number_fragments = 0;//number of fragments to seed in the material

	int m_seed =0;
	//TODO: not currently exposed
	int tries_per_segment = 10;// number of times the algorithm will retry finding an initial segment
	//when all weights of the elemetns are equal 1 will work but this needs increased in porportion to the maximun difference in element weights

	vec3d	vessel_orient_weights = vec3d(0.9090909090, 0.9090909090, 0.0);			// W.x = Weight for collagen orientation
	//W.y = Weight for previous vessel direction
	//W.z is unused
	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 

	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);

	double min_segment_length = 0.1;

	int fragment_seeder = 1;//0 classic fragment seeder, 1 multidomian fragment seeder

	int angio_boundary_type = 1;// 0 same as non angio materials, 1 pass through
	int angio_boundary_groups = 1;//each bit in this parameter fragments can only travel between the groups they are in

	double density_gradient_threshold = 0.01;//set this to be higher to turn off the direction change on encountering different densities
};
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
	virtual bool SeedFragments(SimulationTime& time, Culture * culture)=0;
	FragmentSeeder(CultureParameters * cp, FEAngio & angio) : m_angio(angio), culture_params(cp) {}
	virtual ~FragmentSeeder(){}
protected:
	FEAngio & m_angio;
	CultureParameters * culture_params;
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
	// Seed an initial fragment within the grid
	bool createInitFrag(Segment& Seg, SegGenItem & item, Culture * culture);
	std::vector<FEDomain *> domains;
};

class MDByVolumeFragmentSeeder : public FragmentSeeder
{
public:
	bool SeedFragments(SimulationTime& time, Culture * culture) override;
	MDByVolumeFragmentSeeder(CultureParameters * cp, FEAngio & angio) : FragmentSeeder(cp, angio) {}
private:
	// Seed an initial fragment within the grid
	bool createInitFrag(Segment& Seg, SegGenItem & item, Culture * culture);
	std::vector<FEDomain *> domains;
};

//-----------------------------------------------------------------------------
// The CULTURE class contains all the functions that describe how the 
// SEGMENT class is to grow and orient itself. These functions are the 
// rule-set that arrange the line segments into vascular networks, 
// mimicking angiogenesis. 
class Culture  
{
public:
	Culture(FEAngio& angio, FEAngioMaterial * matl, CultureParameters * cp);
	virtual ~Culture();

	// initialize
	bool Init();

	// Perform a growth step
	void Grow(SimulationTime& time);

	// Reposition the vessels based on the FE solution
	void Update();

	// get the total vessel length
	double TotalVesselLength() const { return m_total_length; }

	// Determine the orientation vector of a newly created segment
	vec3d FindGrowDirection(Segment::TIP& tip) const;
	
	// Find the density-based length scale factor at a point of the grid
	double FindDensityScale(const GridPoint& pt) const;
	
	// Add a segment to the culture (without checking BCs).
	// This assumes that the segment is valid.
	void AddSegment(Segment& seg);

	// Add a new segment to the culture.
	// This will apply BCs to the new segment and may result in 
	// the addition of several new segments. 
	void AddNewSegment(Segment& seg);

	// return the number of segments
	int Segments() const { return m_nsegs; }

	// get the segment list
	const SegmentList& GetSegmentList() const { return m_frag; }

	// return the active segment list
	const SegmentTipList& GetActiveTipList() const { return m_active_tips; }

	// get the total number of active tips
	int ActiveTips() const { return static_cast<int>(m_active_tips.size()); }


	void ChangeBC(FEAngio & angio, int bcset);

	void CreateBranchingForce(Segment& seg);

	// Find the active tips
	void FindActiveTips();

protected:
	virtual void SetWeights(vector<SegGenItem> & weights, std::vector<FEDomain*> & domains);

private:

	// grow vessels
	void GrowVessels();

	// branching phase
	void BranchVessels(SimulationTime& time);

	// create a branch
	void BranchSegment(Segment::TIP& it);

	// fuse segments (i.e. anastomosis)
	void FuseVessels();

	// Grow a segment
	Segment GrowSegment(Segment::TIP& it, bool branch = false, bool bnew_vessel = false, vec3d growthDirection = vec3d());

	// Create a new segment connecting two existing segments that are fusing through anastomosis
	static Segment ConnectSegment(Segment& it, Segment& it2, int k, int kk);

	// Update the new vessel length 
	void UpdateNewVesselLength(SimulationTime& time);

	


	int					m_nsegs;				// Counter that stores in current number of Segments within the simulation domain
	double				m_total_length;			// Total vascular length within the domain (sum of the length of all Segments) (in um)
	list<Segment>		m_frag;					// vessel fragments
	list<Segment::TIP*> m_active_tips;			// list of active tips
	CultureParameters * m_cultParams;
	BC *		bc;
	FEAngio&	m_angio;
	
	FragmentSeeder * fseeder = nullptr;
public:
	double	m_vess_length;	// new segment length
	int		m_num_vessel;   // Counter that indicates the next vessel ID number of new Segments
	int		m_num_branches;		// Counter indicating the number of branches formed during the simulation
	FEAngioMaterial * m_pmat = nullptr;
	
	int m_num_anastom;	// Counter indicating the number of anastomoses formed during the simulation
	int m_num_zdead;

	friend class BouncyBC;
	friend class BC;
	friend class StopBC;
	friend class MBC;
	friend class PassThroughMBC;
};
