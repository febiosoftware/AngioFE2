#pragma once
#include "StdAfx.h"
#include "Fileout.h"
#include "FESproutBodyForce.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngioMaterial.h"
#include "Segment.h"
#include "FECore/FESolidDomain.h" //isd this include correct or should i just forward declare the class



//-----------------------------------------------------------------------------
class FEModel;
class Culture;
class FEAngioMaterial;

//-----------------------------------------------------------------------------
// This class represents the time parameters
class SimulationTime
{
public:
	double	t;		// current time value
	double	dt;		// time increment from last time
	double	maxt;	// end time of simulation
	SimulationTime() { t = dt = maxt = 0.0; }
};

//-----------------------------------------------------------------------------
// The FEAngio class contains all the functionality of the growth model.
class FEAngio
{
public:
	FEAngio(FEModel& fem);
	~FEAngio();

	// initialize the FEAnio class
	bool Init();

	// Get the FE model
	FEModel& GetFEModel() const;

	//check the get FEModel above it may not be useful in any way
	FEMesh * GetMesh() const;

	// get the run time (in secs)
	double RunTime() const;

	// get the simulation time
	SimulationTime& CurrentSimTime() { return m_time; }

	void adjust_mesh_stiffness();
	void update_sprout_stress_scaling();
	static FEAngioMaterial* FindAngioMaterial(FEMaterial* pm);
	double FindECMDensity(const GridPoint& pt);
	vec3d LocalToGlobal(FESolidElement * se, vec3d & rst) const;
	vec3d FindRST(const vec3d & r, vec2d rs, FESolidElement * elem) const;
	GridPoint FindGridPoint(FEDomain * domain, int nelem, vec3d& q) const;
	vec3d Position(const GridPoint& pt)  const;
	// Creates a vector of specified paramenter
	std::vector<double> createVectorOfMaterialParameters(FEElement * elem,
		double FEAngioNodeData::*materialparam);
	// gets the value of a parameter at a given point interpolated from the shape function
	double genericProjectToPoint(FEElement * elem,
		double FEAngioNodeData::*materialparam,const vec3d & pos);
	
	//some funtions to replace the loops everywhere in the code
	//for each excludes non-angio materials
	void ForEachNode(std::function<void(FENode &)> f);
	void ForEachNode(std::function<void(FENode &)> f, std::vector<int> & matls);
	void ForEachElement(std::function<void(FESolidElement&, FESolidDomain&)> f);
	void ForEachElement(std::function<void(FESolidElement&, FESolidDomain&)> f, std::vector<int> & matls);
	void ForEachDomain(std::function<void(FESolidDomain&)> f);
	void ForEachDomain(std::function<void(FESolidDomain&)> f, std::vector<int> & matls);

	vec3d uniformRandomDirection();
	vec3d uniformInUnitCube();

	//calcualtes the gradient at the given natural coordinates
	static vec3d gradient(FESolidElement * se,std::vector<double> & fn, vec3d pt);

	static bool IsInsideHex8(FESolidElement * se, vec3d y, FEMesh * mesh, double r[3]);

	//these freindships are for displaying/reading the data and are okay
	friend class Fileout;
	friend class FEPlotAngioCollagenFibers;
	friend class FEPlotAngioECMDensity;
	friend class FEPlotAngioECMAlpha;
	friend class FEPlotAngioGradient;
	//the following friendships are bad and need removed eventually
	//TODO: remove the freindship, creation in the old way requires this
	//or consider making the node and element data public
	friend class FEAngioMaterial;
	friend class BC;
private:
	
	// Initialize the nodal ECM values
	bool InitECMDensity();
	void UpdateECM();
	//will set the times and whether or not this is the final angio step in this mechanical step
	int FindGrowTimes(std::vector<std::pair<double, double>> & time_pairs, int start_index);

	// Initialize nodal collagen fiber directions
	bool InitCollagenFibers();

	// Initialize nodal collagen fiber directions
	bool InitSoluteConcentration();

	int FindVEGF();

	// Init FE stuff
	bool InitFEM();

	void FinalizeFEM();

	void SetupSurface();

	// do the final output
	void Output();

	static bool feangio_callback(FEModel* pfem, unsigned int nwhen, void* pd)
	{
		FEAngio* pfa = (FEAngio*)(pd);
		pfa->OnCallback(pfem, nwhen);
		return true;
	}

	void OnCallback(FEModel* pfem, unsigned int nwhen);

public:	// parameters read directly from file

	//consider converting these to vectors or arrays as the multidomain simulations require the nde/element ids to be sequential
	std::unordered_map<int, FEAngioNodeData> m_fe_node_data;
	std::unordered_map<int, FEAngioElementData> m_fe_element_data;

	// miscellaneous
	unsigned int	m_irseed;			// Seed number for the random number generator

    double	half_cell_l;			// Half the length of one grid element in the x direction, use in determining variable time step
    
	int		total_bdyf;
	int		FE_state;			// State counter to count the number of solved FE states

	std::vector<FEAngioMaterial*>	m_pmat;	//!< the angio-material pointer
	std::vector<int>                m_pmat_ids;

	angiofe_random_engine rengine;

private:
	FEModel&		m_fem;		// the FE model
	//both nodes and elements id's go from 1 to n+1 for n items
	//first element is padding so the id can be used to lookup the data for that node
	

	SimulationTime	m_time;		// simulation time

    time_t m_start;			// time of start
	Fileout fileout;		// output manager
	
	std::uniform_real_distribution<double> ztopi;
	std::uniform_real_distribution<double> zto2pi;
	std::uniform_real_distribution<double> n1to1;
};
