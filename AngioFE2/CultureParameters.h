#pragma once
#include <FECore/vec3d.h>

const int MAXPARAMSIZE = 256;


//contains the material parameters and the accessor fucntions when the variables are made dependent on time
class CultureParameters
{
public:
	CultureParameters()
	{
		m_boundary_condition_type[0] = 's';
		m_boundary_condition_type[1] = '\0';
		vessel_file[0] = '\0';
	}

	//not currently exposed
	//These parameters govern the amount of growth at each time step for the microvessels
	//          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double	m_y0 = -19.1278;   // Bottom of sigmoid curve
	double	m_culture_a = 1900.0;    // Distance between top and bottom of the curve (a + y0 = top of curve)
	double	m_x0 = 4.9474;   // Time point at which t is halfway between top & bottom
	double	m_culture_b = 1.4549;    // Steepness of the curve

	//this may supercede the above parameters
	double growth_length_over_time = 140.0;//this is the number of units a segment would grow over 1 unit of time if the density is 3mgml. by using a curve the .feb creator has control of the distribution used and its approximation

	//the initial material parameters a,b,N
	double sprout_s_mag = 3.72e-12;
	double sprout_s_range = 0.004;
	double sprout_s_width = 2.0;

	//dont expose
	double	m_initial_vessel_length = 0.0;	// initial vessel length
	double	m_initial_branch_probability = 0.0;     // Probability that initial segments will branch
	bool m_bzfibflat = false;//flatten the fibers in the z direction

	//parameters adjustable by the user for vessel network morphology
	// Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjustment = 1.0;
	// If the shortest distance vector between a segment and any other segment is less than this value
	//TODO: I think making this a percentage of the growth length makes more sense.
	double m_anastomosis_distance = 25;
	bool	m_branch = false; //whether the vessels are allowed to branch
	bool	m_anastomosis = false; //whether the vessels are allowed to fuse together
	// vessel_width - Diameter of microvessels (Default: 7 um)
	double m_vessel_width = 7.0;
	char m_boundary_condition_type[MAXPARAMSIZE];//currently s for stop, or b for bouncy 

	//parameters for ECM density/alignment
	int m_matrix_condition = 0; // flag indicating how the collagen fibers are oriented initially ( 0 = random classic mode, 1 multimaterial mode,  3 = along local element direction)
	int ecm_control = 0; //flag indicating how the ecm density and anisotropy are initialized: 0 constant mode, 1 specified mode, 2 no overwrite 
	double m_matrix_density = 3.0;//mg/ml 3.0 default

	//symmety
	vec3d m_symmetry_plane = vec3d(0.0, 0.0, 0.0); //symmetry plane

	//uncategoriezed
	int m_composite_material = 1;//whether or not the material is a composite material

	//TODO: doesn't appear to be the same as sprout_s_mag
	double m_sprout_force = 1.0;

	int m_number_fragments = 0;//number of fragments to seed in the material

	//TODO: not currently exposed
	int tries_per_segment = 10;// number of times the algorithm will retry finding an initial segment
	//when all weights of the elemetns are equal 1 will work but this needs increased in porportion to the maximun difference in element weights

	vec3d	vessel_orient_weights = vec3d(0.9090909090, 0.9090909090, 0.0);			// W.x = Weight for collagen orientation
	//W.y = Weight for previous vessel direction
	//W.z is unused
	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 

	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);

	double min_segment_length = 0.1;

	int fragment_seeder = 1;//0 classic fragment seeder, 1 multidomian fragment seeder, 2 by volume fragment seeder, 3 from file

	int branching_scheme = 2;

	char vessel_file[MAXPARAMSIZE];//only used if seeding from file

	int angio_boundary_type = 1;// 0 same as non angio materials, 1 pass through
	int angio_boundary_groups = 1;//each bit in this parameter fragments can only travel between the groups they are in


	double density_gradient_threshold = 0.01;//set this to be higher to turn off the direction change on encountering different densities

	double average_length_to_branch_point = 80.0;//sets the average length before a branch is encountered
	double std_deviation = 5.0;
	//consider allowing the user to specify the distribution and any paramters associated with the distribution that is used
	double emerge_time_mean = 0.5;
	double emeerge_time_std_deviation = 0.25;


	const bool io = true;
	friend class FEParamContainer;
	friend class FEAngioMaterial;

	double GetBranchProbility(double dt)const{ return m_branch_chance * dt; }
	double GetWeightInterpolation(double dt) const { return m_weight_interpolation * dt; }
private:
	double	m_branch_chance = 0.1;    // Probability of forming a new branch
	double m_weight_interpolation = 0.4; //used to control the amount of affect of the collagen direction vs the previous direction for an approximation over 1t

};