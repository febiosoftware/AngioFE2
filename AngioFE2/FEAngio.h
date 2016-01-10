///////////////////////////////////////////////////////////////////////
// FEAngio.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The ANGIO class contains all the functionality of the growth model.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "FECore/FEModel.h"
#include "Fileout.h"
#include "FEAngioMaterial.h"
#include "FESproutBodyForce.h"
#include "Grid.h"
#include "Culture.h"

//-----------------------------------------------------------------------------
class FEModel;

//-----------------------------------------------------------------------------
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

public:
	void initBranch();
	void updateTime();
	void updateTotalLength();
	virtual void Growth();
	void updateLength();
	virtual void Branch(list<Segment>::iterator it);
	virtual void Fuse();
    void check4anast(list<Segment>::iterator it, int k);
	//void search_nearby_neighbors(list<Segment>::iterator it, int k, int elem_num);
	//void scan_4_segs(list<Segment>::iterator it, int k, int elem_num);
    void anastomose(double dist0, double dist1, int k, list<Segment>::iterator it, list<Segment>::iterator it2);
	void removeErrors();
	void find_active_tips();
	void kill_dead_segs();

public:
	void Growth(FEModel& fem);
	void Branch(list<Segment>::iterator it, FEModel& fem);
	bool Subgrowth(int sub_steps, FEModel& fem);
	void displace_vessels();
	void apply_sprout_forces(FEModel& fem, int load_curve, double scale);
	int create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, FEModel& fem, int load_curve);
	void update_body_forces(FEModel& fem, double scale);
	void create_branching_force(Segment& seg, FEModel& fem);
	void update_grid(FEMesh &mesh);
	void update_ECM();
	mat3d calculate_deform_tensor(Elem elem, double e1, double e2, double e3);
	vec3d shapefun_d1(const double xix, const double xiy, const double xiz, int node);
	void adjust_mesh_stiffness(FEModel& fem);
	void initialize_grid_volume();
	void update_grid_volume();
	void update_ecm_den_grad();
	void output_params();
	void save_cropped_vessels();
	void update_sprout_stress_scaling();
	void circ_gel();
	void update_angio_sprout(int i, bool bactive, const vec3d& rc, const vec3d& sprout_vect);

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
	double			m_branch_chance;    // Probability of forming a new branch
    double			m_a1,m_a2,m_a3;     // Parameters for the branching curve (TODO: Make these user parameters)

	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 	                                                            
	// TODO: Makes these user parameters
	//          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double m_y0;   // Bottom of sigmoid curve
	double m_a;    // Distance between top and bottom of the curve (a + y0 = top of curve)
	double m_x0;   // Time point at which t is halfway between top & bottom
	double m_b;    // Steepness of the curve
	double m_d;    // Initial value of the curve (t = 0)
    double m_vess_length;		// TODO: Move this to Culture?

    // Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjust;                                       

	
	// If the shortest distance vector between a segment and any other segment is less than this value
	// TODO: I think making this a percentage of the growth length makes more sense.
    double m_anast_dist;

public:
    double	half_cell_l;		// Half the length of one grid element in the x direction, use in determining variable time step
	int		m_sub_cycles;		// number of FE steps per growth step
    int		m_n;				// number of total time steps ?
	double	m_dtA, m_dtB, m_dtC;	// Time step parameters. TODO: make these user variables.
    
	list<list<Segment>::iterator > active_tips;
	
	vector<vector<double> > sprout_nodes;

public:

	int total_bdyf;
	double FE_time_step;
	int FE_state;
	
	bool yes_branching;
	bool yes_anast;

	// time stepping parameters
	double m_maxt;		// Time at which the simulation ends (in days)
	double m_t;			// Current time (in days)
	double m_dt;		// Variable time step (in days)

public:
	// stats
	double	m_total_length;		// Total vascular length within the domain (sum of the length of all Segments) (in um)
	int		m_num_branches;		// Counter indicating the number of branches formed during the simulation
	int		m_nsegs;			// Counter that stores in current number of Segments within the simulation domain
	bool	m_branch;			// Boolean flag that indicates to the model that the Segment being created is the result of a new branch

public:	// IO stuff

    time_t m_start;				// time of start
	Fileout fileout;
	FILE *killed_segs_stream; 
	bool kill_off;
	
public:
	FESproutBodyForce*	m_pbf;	//!< sprout body-force
	FEAngioMaterial*	m_pmat;	//!< the angio-material pointer

private:
	FEModel&		m_fem;	// the FE model
    Grid			m_grid;
	Culture*		m_pCult;  
};
