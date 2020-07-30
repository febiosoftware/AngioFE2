#pragma once
#include <FECore/vec3d.h>

const int MAXPARAMSIZE = 256;


//contains the material parameters and the accessor fucntions when the variables are made dependent on time
class CultureParameters
{
public:
	CultureParameters()
	{
	}

	//this may supercede the above parameters
	double growth_length_over_time = 140.0;//this is the number of units a segment would grow over 1 unit of time if the density is 3mgml. by using a curve the .feb creator has control of the distribution used and its approximation

	//the initial material parameters a,b,N
	double sprout_s_mag = 3.72e-12;
	double sprout_s_range = 0.004;
	double sprout_s_width = 2.0;

	bool m_bzfibflat = false;//flatten the fibers in the z direction

	//parameters adjustable by the user for vessel network morphology
	// Length adjuster, scales the length of new segments
	// allows fine tuning of total vascular length produced by the simulation.
	double m_length_adjustment = 1.0;

	// vessel_width - Diameter of microvessels (Default: 7 um)
	double m_vessel_width = 7.0;

	//parameters for ECM density/alignment
	int ecm_control = 0; //flag indicating how the ecm density and anisotropy are initialized: 0 constant mode, 1 specified mode, 2 no overwrite 
	double m_matrix_density = 3.0;//mg/ml 3.0 default
	int apparent_density = 1; //flag to use only initial density to scale rather than apparent. 0 means use initial density 1 means use apparent (legacy)

	//symmety
	vec3d m_symmetry_plane = vec3d(0.0, 0.0, 0.0); //symmetry plane

	//uncategoriezed
	int m_composite_material = 1;//whether or not the material is a composite material

	//TODO: doesn't appear to be the same as sprout_s_mag
	double m_sprout_force = 1.0;
	double m_sprout_force_x0 = 2.0;
	double m_sprout_force_y0 = -0.004;
	double m_sprout_force_a = 1.0081;
	double m_sprout_force_b = 0.5436;

	//TODO: not currently exposed
	int tries_per_segment = 10;// number of times the algorithm will retry finding an initial segment
	//when all weights of the elemetns are equal 1 will work but this needs increased in porportion to the maximun difference in element weights

	// Microvessel growth rate is modeled using the Boltzmann sigmoid curve 

	vec3d m_density_scale_factor = vec3d(-0.016, 5.1605, 0.5112);

	double min_segment_length = 0.1;

	int active_tip_threshold = 500;
	double stress_radius = 200.0;

	int max_recursion_depth = 5;

	const bool io = true;
	friend class FEParamContainer;
	friend class FEAngioMaterial;

	double GetWeightInterpolation(double dt) const { return m_weight_interpolation * dt; }
private:
	double m_weight_interpolation = 0.4; //used to control the amount of affect of the collagen direction vs the previous direction for an approximation over 1t

};