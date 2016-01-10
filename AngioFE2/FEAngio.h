///////////////////////////////////////////////////////////////////////
// FEAngio.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The ANGIO class contains all the functionality of the growth model.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "FECore/FEModel.h"
#include "Fileout.h"
#include "Profiler.h"
#include "FEAngioMaterial.h"
#include "FESproutBodyForce.h"
#include "Grid.h"
#include "Data.h"
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
	void save_vessel_state();
	void save_bdy_forces(FEModel& fem);
	void save_time();
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

public:
    double	half_cell_l;		// Half the length of one grid element in the x direction, use in determining variable time step
	int		m_sub_cycles;		// number of FE steps per growth step
    
	list<list<Segment>::iterator > active_tips;
	
	vector<vector<double> > sprout_nodes;

    Grid grid;
	Data data;
	Culture cult;  

public:

	int total_bdyf;
	double FE_time_step;
	int FE_state;
	
	bool yes_branching;
	bool yes_anast;

public:	// IO stuff

	Fileout fileout;
	
	FILE *stream;                                                           
	FILE *bf_stream;

	FILE *time_stream;
	bool time_write_headers;

	FILE *killed_segs_stream; 
	bool kill_off;

public:
	FESproutBodyForce*	m_pbf;	//!< sprout body-force
	FEAngioMaterial*	m_pmat;	//!< the angio-material pointer

private:
	FEModel&	m_fem;	// the FE model
};
