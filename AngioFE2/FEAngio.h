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

	// Get the FE model
	FEModel& GetFEModel();

	// Get the culture
	Culture& GetCulture() { return *m_pCult; }

	// Get the grid
	Grid& GetGrid() { return m_grid; }

	// get the run time (in secs)
	double RunTime();

	// get the simulation time
	SimulationTime& CurrentSimTime() { return m_time; }

	// Total number of sprouts
	int Sprouts();

public:
	void adjust_mesh_stiffness();
	void update_ecm_den_grad();
	void update_sprout_stress_scaling();

private:
	// Initialize the nodal ECM values
	bool InitECMDensity();

	// Initialize nodal collagen fiber directions
	bool InitCollagenFibers();

	// Init FE stuff
	bool InitFEM();

	// do the final output
	void Output();

	// create the sprouts
	void CreateSprouts(double scale);

	// update sprouts
	void UpdateSprouts(double scale);

private:
	static void feangio_callback(FEModel* pfem, unsigned int nwhen, void* pd)
	{
		FEAngio* pfa = (FEAngio*)(pd);
		pfa->OnCallback(pfem, nwhen);
	}

	void OnCallback(FEModel* pfem, unsigned int nwhen);

public:	// parameters read directly from file

	// sprout force parameters
	double	m_tip_range;		// sprout force range
	double	m_sproutf;			// sprout force magnitude
	double	m_spfactor;			// Sprout force directional factor
    int		m_bsp_sphere;		// Flag for sprout force representations

	// boundary conditions
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

    double	half_cell_l;			// Half the length of one grid element in the x direction, use in determining variable time step
    int		m_ntime;				// number of total time steps ?
	double	m_dtA, m_dtB, m_dtC;	// Time step parameters. TODO: make these user variables.
    
public:
	int		total_bdyf;
	int		FE_state;			// State counter to count the number of solved FE states

public:
	FESproutBodyForce*	m_pbf;	//!< sprout body-force
	FEAngioMaterial*	m_pmat;	//!< the angio-material pointer

private:
	FEModel&		m_fem;		// the FE model
    Grid			m_grid;		// The grid class
	Culture*		m_pCult;	// The culture class

	SimulationTime	m_time;		// simulation time

    time_t m_start;			// time of start
	Fileout fileout;		// output manager
};
