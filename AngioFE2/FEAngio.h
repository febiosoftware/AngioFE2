#pragma once
#include "Fileout.h"
#include "FESproutBodyForce.h"
#include "Grid.h"

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

public:
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

	// run the simulation
	bool Run();

	// Get the FE model
	FEModel& GetFEModel();

	// Get the culture
	Culture& GetCulture() { return *m_pCult; }

	// Get the grid
	Grid& GetGrid() { return m_grid; }

public:
	void updateTime();
	void removeErrors();

public:
	bool Subgrowth(int sub_steps);
	void apply_sprout_forces(int load_curve, double scale);
	int create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, int load_curve);
	void update_body_forces(double scale);
	void update_ECM();
	void adjust_mesh_stiffness();
	void update_ecm_den_grad();
	void output_params();
	void update_sprout_stress_scaling();
	void circ_gel();
	void update_angio_sprout(int i, bool bactive, const vec3d& rc, const vec3d& sprout_vect);

private:
	// Initialize the nodal ECM values
	bool InitECMDensity();

	// Initialize nodal collagen fiber directions
	bool InitCollagenFibers();

	// Init FE stuff
	bool InitFEM();

	// Run the FE model
	bool RunFEM();

	// Update the grid based on the FE solution.
	void UpdateGrid();

	// Reposition the vessels based on the FE solution
	void DisplaceVessels();

	// Recalculate the total lenght of the network
	void UpdateTotalLength();

	// do the final output (called at the end of Run())
	void Output();

protected:
	// enforce_fiber_BCS - If a node lies on a boundary face, adjust the collagen fiber orientation at that node accordingly
	void enforce_fiber_BCS(Node &node, bool circ);

public:	// parameters read directly from file

	// sprout force parameters
	double	m_tip_range;		// sprout force range
	double	m_sproutf;			// sprout force magnitude
	double	m_spfactor;			// Sprout force directional factor
    int		m_bsp_sphere;		// Flag for sprout force representations

	// boundary conditions
	char m_cgelbc;					// Boundary conditions for the gel (LAC, SAC, UNC)
	double m_Sx;						// Location of the symmetry plane along the x-direction
	double m_Sy;						// Location of the symmetry plane along the y-direction
	double m_Sz;						// Location of the symmetry plane along the z-direction

	// miscellaneous
	unsigned int	m_irseed;			// Seed number for the random number generator
	int				comp_mat;			// is composite material used (TODO: do we still need this?)
	double			phi_stiff_factor;	// stiffness factor, which scales the amount of mesh displacement that's sent to the microvessels
	int				m_bsprout_verify;	// Flag for the sprout verification problem
    double			m_vessel_width;     // Width of Segments (in um) (TODO: make this a user parameter)
	int				m_matrix_cond;		// flag indicating how the collagen fibers are oriented initially ( 0 = random, 3 = along local element direction)
	int				m_bzfibflat;		// flatten fiber option

public:
    double	half_cell_l;			// Half the length of one grid element in the x direction, use in determining variable time step
	int		m_sub_cycles;			// number of FE steps per growth step
    int		m_ntime;				// number of total time steps ?
	double	m_dtA, m_dtB, m_dtC;	// Time step parameters. TODO: make these user variables.
    
	vector<vector<double> > sprout_nodes;

public:

	int		total_bdyf;
	int		FE_state;			// State counter to count the number of solved FE states
	
	// time stepping parameters
	SimulationTime	m_time;

public:
	// stats
	double	m_total_length;		// Total vascular length within the domain (sum of the length of all Segments) (in um)

public:	// IO stuff

    time_t m_start;				// time of start
	Fileout fileout;
	bool kill_off;
	
public:
	FESproutBodyForce*	m_pbf;	//!< sprout body-force
	FEAngioMaterial*	m_pmat;	//!< the angio-material pointer

private:
	FEModel&		m_fem;		// the FE model
    Grid			m_grid;		// The grid class
	Culture*		m_pCult;	// The culture class
};
