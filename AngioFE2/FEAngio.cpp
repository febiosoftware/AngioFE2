///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include <FECore/log.h>
#include "FECore\FESolidDomain.h"
#include "FEBioMech\FEElasticMaterial.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "Elem.h"
#include "angio3d.h"
#include <iostream>

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, FEMaterial* pmat);

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : m_fem(fem)
{
	cout << endl << "Angiogenesis Growth Model:" << endl << endl;

	yes_branching = true;														// Flag to turn branching on/off
	yes_anast = true;															// Flag to turn anastomosis on/off
	
	total_bdyf = 0;													// Body force counter
	
	FE_state = 0;													// State counter to count the number of solved FE states
	FE_time_step = 0.5;												// Time step within FEBio

	comp_mat = 0;
	m_pmat = 0;

	phi_stiff_factor = 1.0;
	m_sub_cycles = 2;
    m_n = 1;
	m_branch = false;

	m_length_adjust = 1.0;
	m_anast_dist = 75.0;

	// initialize time stepping parameters
	m_t = 0.0;
    m_maxt = 0.0;
	m_dt = 0.25;

	// TODO: What are these and make these user parameters
    m_dtA = 0.0637;
    m_dtB = 9.0957;
    m_dtC = 2.6073;


	// a1, a2, a3 - Parameters for branching curve
    m_a1 = -1.2653;
	m_a2 = 1.535;
	m_a3 = 0.1108;

	// parameter for growth curve (TODO: Maybe move all this to Culture?)
    m_a = 1900.0;
    m_b = 1.4549;
    m_x0 = 4.9474;
    m_y0 = -19.1278;
    m_d = m_y0 + m_a/(1+pow(E,m_x0/m_b));                                 // d - Initial value of growth curve (t = 0)

	m_vess_length = m_d;

	// vessel_width - Diameter of microvessels (Default: 7 um)
	m_vessel_width = 7;

	// branch probability
	m_branch_chance = 0.1;

	m_nsegs = 0;			// nsegs - Initialize segment counter
    m_total_length = 0.;
	m_num_branches = 0;		// Initialize branching counter
    
	m_bsprout_verify = 0;				// Sprout verification problem flag

	// boundary conditions
	m_cgelbc = 'u';					// Gel boundary conditions ('u' unconstrained)
	m_Sx = 0.;						// Location of the x-symmetry plane
	m_Sy = 0.;						// Location of the y-symmetry plane
	m_Sz = 0.;						// Location of the z-symmetry plane

	// Input random seed number
	m_irseed = (unsigned int) time(0);

	// default sprout force parameters
   	m_sproutf   = 1.0;		// Sprout force magnitude	
	m_tip_range = 250.;		// Sprout force range
	m_spfactor = 0.;		// Sprout force directional factor
	m_bsp_sphere = 0;		// Switch between local directional (0), local isotropic (1), and global isotropic (2) sprout froce representations

	kill_off = false;
	killed_segs_stream = fopen("out_dead_segs.ang","wt");
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
	fclose(killed_segs_stream);
}

//-----------------------------------------------------------------------------
FEModel& FEAngio::GetFEModel()
{
	return m_fem;
}

//-----------------------------------------------------------------------------
// find the angio material component
FEAngioMaterial* FindAngioMaterial(FEMaterial* pm)
{
	FEMaterial* pmat = pm->FindComponentByType("angio");
	if (pmat)
	{
		FEAngioMaterial* pma = dynamic_cast<FEAngioMaterial*>(pmat);
		return pma;
	}
	return 0;
}

//-----------------------------------------------------------------------------
// Initializes the FEAngio object.
bool FEAngio::Init()
{
	// Seed the random number generator
	srand(m_irseed);

	// Print out the seed number for the random generator
	fileout.printrandseed(m_irseed);					

	// If branching is turned off, we set the branching chance to zero
	if (yes_branching == false) m_branch_chance = 0.0;

	vector<vec3d> fiber;
	if (grid.load_cond == 3)
	{
		// get the material
		FEMaterial* pm = m_fem.GetMaterial(0);
		FEMaterial* efd = pm->FindComponentByType("EFD neo-Hookean");
		if (efd == 0) grid.load_cond = 0;
		else
		{
			if (CreateFiberMap(fiber, efd) == false) grid.load_cond = 0;
		}
	}

	vector<double> density;
	if(grid.coll_den == 0.0)
	{
		// get the material
		FEMaterial* pm = m_fem.GetMaterial(0);
		FEMaterial* am = pm->FindComponentByType("angio");

		if (am == 0) grid.coll_den = 3.0;
		else
		{
			if (CreateDensityMap(density, am) == false) grid.coll_den = 3.0;
		}

	}

	// create the grid based on the FEBio mesh
	grid.CreateGrid(m_fem.GetMesh(), fiber, density);

	// start timer
	time(&m_start);
	
	//// FEBIO - Set parameters for the finite element model
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	pstep->m_dt = pstep->m_dt0;							// Set the current time step to the initial time step
	pstep->m_ntimesteps = 0;									// Set the number of time steps to 0
	pstep->m_iteopt = 100;									// Set the optimal iterations to 100
	pstep->m_maxretries = 10;									// Set the maximum retries for the autotimestepper to 10
	pstep->m_dt = FE_time_step;

	// See if an "angio" material is defined.
	bool bmat = true;
	FEAngioMaterial* pma = FindAngioMaterial(m_fem.GetMaterial(0));
	if (pma == 0) bmat = false;
	else 
	{
		//pma->scale = (1.0/4.0)*0.001;

		if (m_Sx != 0.){ pma->Sx = m_Sx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (m_Sy != 0.){ pma->Sy = m_Sy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (m_Sz != 0.){ pma->Sz = m_Sz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pma->ApplySym();

		m_pmat = pma;
		update_sprout_stress_scaling();
	}

	// If the angio material is not defined we apply the "old" body force approach
	if (bmat == false)
	{
		FESproutBodyForce* pbf = new FESproutBodyForce(&m_fem);			// Define the one-and-only bodyforce	
		m_pbf = pbf; // --- SAM ---
	
		m_fem.AddBodyLoad(pbf);										// Add the body force to the FEmodel
		FEParameterList& pl = pbf->GetParameterList();					// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);							// Get the magnitude parameter
		FEParam* pb = pl.Find("b"); assert(pb);							// Get the range parameter
	
		pa->value<double>() = (1.0/4.0)*0.001*m_sproutf;	// Set the mangnitude parameter using the input file
		pb->value<double>() = 1.0/m_tip_range;				// Set the range parameter using the input file
		pbf->m_factor = m_spfactor;							// Set the directional factor parameter using the input file
	
		if (m_bsp_sphere == 1) pbf->m_factor = 0.;				// If using local isotropic sprout force, set directional factor to 0

		if (m_bsp_sphere == 2){									// If using global isotropic sprout froce, set directional factor and range to 0 
			pbf->m_factor = 0.; 
			pb->value<double>() = 0;}					

		//cout << pbf->m_factor << endl;									// Print out the directional factor parameter

		if (m_Sx != 0.){ pbf->Sx = m_Sx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (m_Sy != 0.){ pbf->Sy = m_Sy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (m_Sz != 0.){ pbf->Sz = m_Sz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pbf->ApplySym();												// Apply any symmetry to the bodyforce class
	}

//	if (m_bcirc)
//		circ_gel();

	//// ANGIO3D - Seed initial vessel frags
	initialize_grid_volume();								// Initialize the volume of each element in the grid
	cult.SeedFragments();							// Create initial microvessel fragments                               
	kill_dead_segs();										// Remove buggy segments
	find_active_tips();										// Update the active tip container

	
	//// FEBIO - Initialize the FE model
	m_pmat->SetGrid(&grid);
	FE_state++;												// Update the FE state

	// report if the stress or body force approach will be used
	if (bmat)
		felog.printf("Angio material found. Stress approach will be used.");
	else
		felog.printf("Angio materia NOT found. Body-force appraoch will be used.");

	return true;
}


///////////////////////////////////////////////////////////////////////
// Member Functions:
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// FEAngio - Growth
//      Vessel elongation is represented by the addition of new line segments at the locations of the active sprouts.
///////////////////////////////////////////////////////////////////////

void FEAngio::Growth(FEModel& fem)
{
    int k;															// Iterator for the segment tips
    list<Segment>::iterator it;										// Declare iterator through list container 'it'
    	
	//// Elongate the active vessel tips
	for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)					// Iterate through each vessel segment...
	{
		for (k=0; k<=1; ++k)											// Iterate through each tip of the vessel segment...
		{
			if (it->m_tip[k].active != 0)											// If tip is active (not 0)...
			{
		        Segment seg;													// Declare SEGMENT object 'seg'
				seg = cult.createNewSeg(it,k);					// Create new vessel segment at the current tip existing segment 

				cult.m_frag.push_front (seg);											// Append new segment at the top of the list 'frag'
				
				++m_nsegs;													// Iterate the total segment counter
			
			}
            
	    }
		    
		if (yes_branching == true)										// If branching is turned on...
			Branch(it, fem);												// Determine if the segment forms a new branch
		    
    }

	if (yes_anast == true)											// If anatomosis is turned on...
		Fuse();															// Determine which segments form anatomoses

	   
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - Branch
//		Branching is modeled as a random process, determine is the segment passed to this function forms a branch
///////////////////////////////////////////////////////////////////////

void FEAngio::Branch(list<Segment>::iterator it, FEModel& fem)
{
    int k;
    //list<Segment>::iterator it;                                         // Declare iterator through list container 'it'
    
    // Generate a random number between 0 and 1. If that number is less than the branching probability at
    // the current time, or initial branching is activated (data.ini_branch = true), then branching occurs
	
    double den_scale = 1.0;												// Declare the density scaling factor

	double xpt = (it->m_tip[1].rt.x + it->m_tip[0].rt.x)/2;								// Set the branch point to be half-way along the length of the segment
	double ypt = (it->m_tip[1].rt.y + it->m_tip[0].rt.y)/2;
	double zpt = (it->m_tip[1].rt.z + it->m_tip[0].rt.z)/2;
	
	den_scale = cult.findDenScale(xpt, ypt, zpt);					// Set the density scaling factor
	

	if ( float(rand())/RAND_MAX < den_scale*m_dt*m_branch_chance/m_t || (it->init_branch == true) )	// The segment generates a random number, if this number is less than the branching probability, or if the segment has an initial branch, the form a branch
    {
	    if ((it->BCdead == 0) && (it->anast == 0))                  // Segments that have encountered a boundary condition or formed an anastomoses may not form a branch
		{                                                           
	        
			Segment seg;                                            // Declare SEGMENT object 'seg'
			m_num_branches = m_num_branches + 1;              // Iterate the total number of branches +1
			m_branch = true;                                     // Branching flag set to 'true.' This tells the program that the
			                                                        // new vessel segment being created is arising from a branch      
    		it->init_branch = false;								// Turn off the initial branch flag
    				
			if (float(rand())/RAND_MAX < 0.5)                       // Randomly determine which node of the parent segment forms the branch
			    k = 0;
		    else
			    k = 1;
    							
			it->m_tip[k].active = sign(0.5f - float(rand())/RAND_MAX);        // Randomly assign the branch to grow as +1 or -1
    				
			seg = cult.createNewSeg(it,k);           // Create the new vessel segment

            cult.m_num_vessel = cult.m_num_vessel + 1;					// Iterate the vessel counter
		    seg.vessel = cult.m_num_vessel;							// Set the segments ID number
    				                    
			it->Recent_branch = 1;
			seg.Recent_branch = 1;                                  // Indicate that this segment just branched, prevents future segments from branching again too soon
    		
			create_branching_force(seg, fem);						// Create a new sprout force for the branch

			cult.m_frag.push_front (seg);                                  // Append new segment at the top of the list 'frag'
			m_branch = false;                                    // Turn off branching flag once branching algorithm is complete
        }
	}                                                              
	

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - Subgrowth
//		Step through growth through a series of sub-steps to produce smoother model results
///////////////////////////////////////////////////////////////////////

bool FEAngio::Subgrowth(int sub_steps, FEModel& fem)
{
	double subgrowth_scale = 0.;									// The subgrowth scaling factor					
	//double mag;
	list<list<Segment>::iterator >::iterator tip_it;				// Iterator through the list of active segment tips
	list<Segment>::iterator frag_it;								// Iterator through the segment container

	FEMesh& mesh = fem.GetMesh();										// Obtain the FE mesh

	update_sprout_stress_scaling();

	for (int k = 1; k <= sub_steps; k++)							// Iterate through the number of substeps...
	{
		subgrowth_scale = ((double)k/(double)sub_steps);				// Update the value of the subgrowth scaling factor
		
		for (tip_it = active_tips.begin(); tip_it != active_tips.end(); ++tip_it)	// For each of the active growth tips...
		{
			Segment& seg = (*(*tip_it));											// Dereference the tip iterator to obtain the active segment

			// Step growth for the active segment
			if (seg.m_tip[0].active == -1) seg.m_tip[0].rt = seg.m_tip[1].rt + seg.uvect*(subgrowth_scale*seg.length);
			if (seg.m_tip[1].active ==  1) seg.m_tip[1].rt = seg.m_tip[0].rt + seg.uvect*(subgrowth_scale*seg.length);
		}

		update_body_forces(fem, 1.0);									// Update the positions of the body forces

		FE_state++;														// Increase the FE state

		// Reset some parameters for FEBio
		FEAnalysis* pstep = fem.GetCurrentStep();
		pstep->m_dt = pstep->m_dt0;
		pstep->m_ntimesteps = 0;
		pstep->m_iteopt = 100;
		pstep->m_maxretries = 10;

		// Solve the FE problem
		if (fem.Solve() == false) return false;

		update_grid(mesh);												// Update the growth model mesh using the FE solution
		displace_vessels();												// Update microvessel position and orientation using the displacement field from the FE solution
		
		fileout.save_vessel_state(*this);											// Save the current vessel state
	}

	
	for (frag_it = cult.m_frag.begin(); frag_it != cult.m_frag.end(); ++frag_it)      // Iterate through all segments in frag list container (it)                               
	{
		frag_it->findlength();												// Calculate the new length of the microvessel segments										
	}
	
	removeErrors();														// Remove bad segments

	return true;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - displace_vessels
//		Use the displacement field from the FE solution to update microvessels into the current configuration
///////////////////////////////////////////////////////////////////////

void FEAngio::displace_vessels()
{
	int k = 0;													// Iterator for segment tips
	int elem_num = 0;											// Element number
	list<Segment>::iterator it;									// Iterator for the segment container FRAG
	double xpt = 0.; double ypt = 0.; double zpt = 0.;			// Position in Cartesian coordinates
	double xix = 0.; double xiy = 0.; double xiz = 0.;			// Position in the element's natural coordinates		
	double shapeF[8] = {0.};									// Array containing the shape function values at the segment's position
	vec3d disp; vec3d weighted_disp;							// Displacement and weighted displacement vectors						
	Elem elem;													// Element containing the segment
    
	for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		for (k=0; k<=1;++k)                                         // Iterate through both segment tips (k)
		{
			xpt = it->m_tip[k].rt.x;											// Set position to the current segment tip
			ypt = it->m_tip[k].rt.y;
			zpt = it->m_tip[k].rt.z;
		    
			int BC_face = 0;
			elem_num = it->m_tip[k].elem;								// Find the element that contains the segment tip
		    
			if (elem_num != -1)										// If the segment is inside the mesh and has a real element number...
			{
				elem = grid.ebin[elem_num];								// Obtain the element
		    
				grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);	// Convert the position to natural coordinates
		    
				grid.shapefunctions(shapeF, xix, xiy, xiz);				// Obtain the values of the shape functions at this position
		    
				// Calculate the displacement vector by interpolating nodal displacement to the segment tip
				disp.x = shapeF[0]*(*elem.n1).u.x + shapeF[1]*(*elem.n2).u.x + shapeF[2]*(*elem.n3).u.x + shapeF[3]*(*elem.n4).u.x + shapeF[4]*(*elem.n5).u.x + shapeF[5]*(*elem.n6).u.x + shapeF[6]*(*elem.n7).u.x + shapeF[7]*(*elem.n8).u.x;
				disp.y = shapeF[0]*(*elem.n1).u.y + shapeF[1]*(*elem.n2).u.y + shapeF[2]*(*elem.n3).u.y + shapeF[3]*(*elem.n4).u.y + shapeF[4]*(*elem.n5).u.y + shapeF[5]*(*elem.n6).u.y + shapeF[6]*(*elem.n7).u.y + shapeF[7]*(*elem.n8).u.y;
				disp.z = shapeF[0]*(*elem.n1).u.z + shapeF[1]*(*elem.n2).u.z + shapeF[2]*(*elem.n3).u.z + shapeF[3]*(*elem.n4).u.z + shapeF[4]*(*elem.n5).u.z + shapeF[5]*(*elem.n6).u.z + shapeF[6]*(*elem.n7).u.z + shapeF[7]*(*elem.n8).u.z;
		    		    
				weighted_disp = disp*phi_stiff_factor;					// Calculate the weighted displacement vector
		
				// store the old position
				vec3d rt_old = it->m_tip[k].rt;
				
				// Update the segment tip position using the new weighted displacement vector
				it->m_tip[k].rt = rt_old + weighted_disp;
			
				if (grid.findelem(it->m_tip[k].rt.x,it->m_tip[k].rt.y,it->m_tip[k].rt.z) == -1){	// If using the weighted displacement vector causes the segment to move outside the mesh...
					weighted_disp = disp;									// Update using the full displacement vector instead
					it->m_tip[k].rt = rt_old + weighted_disp;
				}
			}
			else													// If the segment doesn't have a real element number...
			{
				it->mark_of_death = true;								// Kill the segment
				it->death_label = 9;
			}			
		}
		
		it->findunit();											// Recalculate the segment's unit vector based on it's new position
		//it->findphi();
	}

    return;
}   


///////////////////////////////////////////////////////////////////////
// FEAngio - apply_sprout_forces
//		Apply sprout forces to the mesh for each active vessel tip
///////////////////////////////////////////////////////////////////////

void FEAngio::apply_sprout_forces(FEModel& fem, int load_curve, double scale)
{
	list<list<Segment>::iterator >::iterator tip_it;				// Iterator for the container that stores active growth tips
	Segment seg;													// Segment placeholder
	vec3d sprout_vect;												// Sprout vector

	vec3d tip(0,0,0);
	double magnitude = scale*m_sproutf;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_t < 4.0)
		magnitude = (1.0/4.0)*m_t*scale;

	//#pragma omp parallel for
	for (tip_it = active_tips.begin(); tip_it != active_tips.end(); ++tip_it)		// For each active growth tip...
	{
		seg = (*(*tip_it));												// Obtain the growth tip

		if (seg.m_tip[0].active == -1){											// If it's a -1 tip...
			tip = seg.m_tip[0].rt;												// Obtain the position of the active tip
			
			// TODO: Why is this calculated? This is already stored in Segment
			sprout_vect = seg.m_tip[0].rt - seg.m_tip[1].rt;							// Calculate the directional unit vector of the sprout
			sprout_vect = sprout_vect/sprout_vect.norm();

			(*tip_it)->m_tip[0].bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, fem, load_curve);}				// Create a new body force, set the tips body force ID
		
		if (seg.m_tip[1].active == 1){											// If it's a +1 tip...
			tip = seg.m_tip[1].rt;												// Obtain the position of the active tip
			
			sprout_vect.x = seg.m_tip[1].rt.x - seg.m_tip[0].rt.x;							// Calculate the directional unit vector of the sprout
			sprout_vect.y = seg.m_tip[1].rt.y - seg.m_tip[0].rt.y;
			sprout_vect.z = seg.m_tip[1].rt.z - seg.m_tip[0].rt.z;
			sprout_vect = sprout_vect/sprout_vect.norm();

			(*tip_it)->m_tip[1].bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, fem, load_curve);}				// Create a new body force, set the tips body force ID
	}

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - create_body_force
//		Add a new body force entry into the body force field applyied to the mesh
///////////////////////////////////////////////////////////////////////

int FEAngio::create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, FEModel& fem, int load_curve)
{
	total_bdyf++;							// Iterate the total body force counter							


//--> SAM
	if (m_pmat)
	{
		m_pmat->AddSprout(vec3d(xpt, ypt, zpt), sprout_vect);
		return m_pmat->Sprouts() - 1;
	}
	else
	{
//<-- SAM
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(fem.GetBodyLoad(0));					// Obtain the body force class
		pbf->AddSprout(vec3d(xpt, ypt, zpt),sprout_vect);												// Add a new component to the body force for this active sprout tip
		return pbf->Sprouts() - 1;																		// Return the ID number for the body force
	}
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_body_forces
//		Update the sprout forces after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_body_forces(FEModel& fem, double scale)
{
	list<Segment>::iterator frag_it;								// Iterator for the segment container FRAG
	vec3d sprout_vect;												// Sprout direction vector

	vec3d tip(0,0,0);
	double magnitude = scale*m_sproutf;								// Magnitude of the sprout force

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_t < 4.0)
		magnitude = (1.0/4.0)*m_t*scale;

	if (m_pmat)
	{
		//m_pmat->scale = magnitude;
		int NSP = m_pmat->Sprouts();
		for (int i=0; i<NSP; ++i)
		{
			FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(i);
			sp.bactive = false;
		}
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(fem.GetBodyLoad(0));			// Obtain the sprout body force field
		int NSP = pbf->Sprouts();										// Obtain the number of sprouts
		for (int i = 0; i < NSP; i++)									// Deactivate all sprout force components
		{
			FESproutBodyForce::SPROUT& sp = pbf->GetSprout(i);				// Obtain the sprout force component
			sp.active = false;												// Deactive the sprout force component
		}
		FEParameterList& pl = pbf->GetParameterList();										// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);												// Get the sprout force magnitude parameter
		pa->value<double>() = magnitude*m_sproutf;													// Set the sprout force magnitude parameter
	}

	FEMesh& mesh = fem.GetMesh();									// Obtain the FE mesh

	//#pragma omp parallel for
	for (frag_it = cult.m_frag.begin(); frag_it != cult.m_frag.end(); ++frag_it)		// Iterate through each segment in the model...
	{
		const Segment& seg = (*frag_it);								// Obtain the segment, keep it constant to prevent changes

		if (((seg.m_tip[0].active == -1) || (seg.m_tip[0].BC == 1)) && (seg.m_tip[0].bdyf_id >= 0)){		  // Turn on the body force for any active -1 segment OR -1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[0] == -1) && (seg.bdyf_id[0] >= 0)){									// Turn on the body force for any active -1 segment
			tip = seg.m_tip[0].rt;																	// Obtain the tip position

			sprout_vect = seg.m_tip[0].rt - seg.m_tip[1].rt;												// Calculate the sprout directional vector
			sprout_vect = sprout_vect/sprout_vect.norm();			

//--> SAM
			update_angio_sprout(seg.m_tip[0].bdyf_id, true, tip, sprout_vect);
//<-- SAM
			}
		
		if (((seg.m_tip[1].active == 1) || (seg.m_tip[1].BC == 1)) && (seg.m_tip[1].bdyf_id >= 0)){		  // Turn on the body force for any active +1 segment OR +1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[1] == 1) && (seg.bdyf_id[1] >= 0)){									// Turn on the body force for any active +1 segment
			tip = seg.m_tip[1].rt;																	// Obtain the tip position
			
			sprout_vect = seg.m_tip[1].rt - seg.m_tip[0].rt;												// Calculate the sprout directional vector
			sprout_vect = sprout_vect/sprout_vect.norm();

//--> SAM
			update_angio_sprout(seg.m_tip[1].bdyf_id, true, tip, sprout_vect);
//<-- SAM
			}
	}
	
	return;
}

//--> SAM
void FEAngio::update_angio_sprout(int id, bool bactive, const vec3d& rc, const vec3d& sprout_vect)
{
	if (m_pmat)
	{
		FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(id);
		sp.bactive = true;
		sp.sprout = sprout_vect;
		m_pmat->UpdateSprout(sp, rc);
	}
	else
	{
		FESproutBodyForce::SPROUT& sp = m_pbf->GetSprout(id);		// Obtain the sprout component 
		sp.active = true;											// Set the sprout component to active
		sp.rc = rc;													// Set the tip position
		sp.sprout = sprout_vect;									// Set the sprout force directional vector
	}
}
//<-- SAM

///////////////////////////////////////////////////////////////////////
// FEAngio - create_branching_force
//		Create a new sprout force component for a newly formed branch		
///////////////////////////////////////////////////////////////////////

void FEAngio::create_branching_force(Segment& seg, FEModel& fem)
{
	vec3d sprout_vect;															// Sprout for directional vector

	total_bdyf = 0;																// Obtain the total number of sprouts
	if (m_pmat) total_bdyf = m_pmat->Sprouts();
	else total_bdyf = m_pbf->Sprouts();
	
	vec3d tip(0,0,0);
	if (seg.m_tip[0].active == -1){														// If the new branch is a -1 segment... 														
		tip = seg.m_tip[0].rt;															// Obtain the position of the new tip
		
		sprout_vect = seg.m_tip[0].rt - seg.m_tip[1].rt;										// Calculate the sprout directional unit vector									
		sprout_vect = sprout_vect/sprout_vect.norm();
		
		seg.m_tip[0].bdyf_id = total_bdyf - 1;}											// Assign the body force ID
		
	if (seg.m_tip[1].active == 1){														// If the new branch is a +1 segment...
		tip = seg.m_tip[1].rt;															// Obtain the position of the new tip
		
		sprout_vect = seg.m_tip[1].rt - seg.m_tip[0].rt;										// Calculate the sprout directional unit vector									
		sprout_vect = sprout_vect/sprout_vect.norm();
		
		seg.m_tip[1].bdyf_id = total_bdyf - 1;}											// Assign the body force ID

//--> SAM
	if (m_pmat)
	{
		m_pmat->AddSprout(tip, sprout_vect);
		total_bdyf = m_pmat->Sprouts();
	}
	else
	{
//!<-- SAM
		m_pbf->AddSprout(tip, sprout_vect);						// Add the new sprout component to the sprout force field
		total_bdyf = m_pbf->Sprouts();												// Update the total number of sprouts
	}
	
	return;
}

//-----------------------------------------------------------------------------
// Update the grid after a deformation using the FE mesh
// TODO: move to the Grid class
void FEAngio::update_grid(FEMesh &mesh)
{
	// loop over all nodes
	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		// Calculate the displacement vector by finding the difference between the node in the deformed FE mesh and the undeformed grid
		grid.nodes[i].u = mesh.Node(i).m_rt - grid.nodes[i].rt;		
				
		// Update the grid node to the current position of the FE mesh
		grid.nodes[i].rt = mesh.Node(i).m_rt;
	}
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_ECM
//		Update the ECM field after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_ECM()
{
	mat3d F;															// Deformation gradient tensor				
	double Jacob = 0.;												// Jacobian (i.e., determinant of F)
	
	vec3d coll_fib;													// Collagen fiber vector
	double ecm_den = 0.;											// Collagen density

	double e1 = 0.; double e2 = 0.; double e3 = 0.;					// Natural coordinates of each node within the element

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i)								// For each element within the mesh...
	{
		Elem& elem = grid.ebin[i];											// Obtain the element
		
		for (int j = 1; j < 9; j++)										// For each node in the element...
		{
			if (j == 1){													// Lower front left node
				ecm_den = (*elem.n1).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n1).collfib.x,(*elem.n1).collfib.y,(*elem.n1).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = -1.;}
			
			if (j == 2){													// Lower front right node
				ecm_den = (*elem.n2).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n2).collfib.x,(*elem.n2).collfib.y,(*elem.n2).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = -1.;}

			if (j == 3){													// Lower back left node
				ecm_den = (*elem.n3).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n3).collfib.x,(*elem.n3).collfib.y,(*elem.n3).collfib.z);	// Set the collagen fiber vector	
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = -1.;}

			if (j == 4){													// Lower back right node
				ecm_den = (*elem.n4).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n4).collfib.x,(*elem.n4).collfib.y,(*elem.n4).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = -1.;}

			if (j == 5){													// Upper front left node
				ecm_den = (*elem.n5).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n5).collfib.x,(*elem.n5).collfib.y,(*elem.n5).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = 1.;}

			if (j == 6){													// Upper front right node
				ecm_den = (*elem.n6).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n6).collfib.x,(*elem.n6).collfib.y,(*elem.n6).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = 1.;}

			if (j == 7){													// Upper back left node
				ecm_den = (*elem.n7).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n7).collfib.x,(*elem.n7).collfib.y,(*elem.n7).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = 1.;}

			if (j == 8){													// Upper back right node
				ecm_den = (*elem.n8).ecm_den0;									// Set the matrix density
				coll_fib = vec3d((*elem.n8).collfib.x,(*elem.n8).collfib.y,(*elem.n8).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = 1.;}
			
			F = calculate_deform_tensor(elem, e1, e2, e3);					// Calculate the deformation gradient tensor
			Jacob = F.det();												// Calculate the Jacobian by taking the determinant of F
			
			coll_fib = F*coll_fib;											// Update the collagen fiber orientation vector into the current configuration using F		

			if (coll_fib.norm() != 0.)										// Normalize the collagen fiber vector to obtain the unit vector
				coll_fib = coll_fib/coll_fib.norm();						

			if (Jacob != 0.)												// Update matrix density using the Jacobian
				ecm_den = ecm_den/Jacob;		
			
			if (j == 1){													// Lower front left node
				if ((*elem.n1).updated == false){								// If the node hasn't been updated yet...
					(*elem.n1).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n1).collfib.y = coll_fib.y;
					(*elem.n1).collfib.z = coll_fib.z;
					(*elem.n1).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n1).updated = true;}										// Set the updated flag
				else if ((*elem.n1).updated == true){							// If the node has been updated...
					(*elem.n1).collfib.x = ((*elem.n1).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n1).collfib.y = ((*elem.n1).collfib.y + coll_fib.y)/2;
					(*elem.n1).collfib.z = ((*elem.n1).collfib.z + coll_fib.z)/2;				
					(*elem.n1).ecm_den = ((*elem.n1).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 2){													// Lower front right node
				if ((*elem.n2).updated == false){								// If the node hasn't been updated yet...
					(*elem.n2).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n2).collfib.y = coll_fib.y;
					(*elem.n2).collfib.z = coll_fib.z;					
					(*elem.n2).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n2).updated = true;}										// Set the updated flag
				else if ((*elem.n2).updated == true){							// If the node has been updated...
					(*elem.n2).collfib.x = ((*elem.n2).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n2).collfib.y = ((*elem.n2).collfib.y + coll_fib.y)/2;
					(*elem.n2).collfib.z = ((*elem.n2).collfib.z + coll_fib.z)/2;
					(*elem.n2).ecm_den = ((*elem.n2).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 3){													// Lower back left node
				if ((*elem.n3).updated == false){								// If the node hasn't been updated yet...
					(*elem.n3).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n3).collfib.y = coll_fib.y;
					(*elem.n3).collfib.z = coll_fib.z;						
					(*elem.n3).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n3).updated = true;}										// Set the updated flag
				else if ((*elem.n3).updated == true){							// If the node has been updated...
					(*elem.n3).collfib.x = ((*elem.n3).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n3).collfib.y = ((*elem.n3).collfib.y + coll_fib.y)/2;
					(*elem.n3).collfib.z = ((*elem.n3).collfib.z + coll_fib.z)/2;
					(*elem.n3).ecm_den = ((*elem.n3).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 4){													// Lower back right node
				if ((*elem.n4).updated == false){								// If the node hasn't been updated yet...
					(*elem.n4).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n4).collfib.y = coll_fib.y;
					(*elem.n4).collfib.z = coll_fib.z;						
					(*elem.n4).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n4).updated = true;}										// Set the updated flag
				else if ((*elem.n4).updated == true){							// If the node has been updated...
					(*elem.n4).collfib.x = ((*elem.n4).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n4).collfib.y = ((*elem.n4).collfib.y + coll_fib.y)/2;
					(*elem.n4).collfib.z = ((*elem.n4).collfib.z + coll_fib.z)/2;
					(*elem.n4).ecm_den = ((*elem.n4).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 5){													// Upper front left node
				if ((*elem.n5).updated == false){								// If the node hasn't been updated yet...
					(*elem.n5).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n5).collfib.y = coll_fib.y;
					(*elem.n5).collfib.z = coll_fib.z;
					(*elem.n5).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n5).updated = true;}										// Set the updated flag
				else if ((*elem.n5).updated == true){							// If the node has been updated...
					(*elem.n5).collfib.x = ((*elem.n5).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n5).collfib.y = ((*elem.n5).collfib.y + coll_fib.y)/2;
					(*elem.n5).collfib.z = ((*elem.n5).collfib.z + coll_fib.z)/2;
					(*elem.n5).ecm_den = ((*elem.n5).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
		
			if (j == 6){													// Upper front right node
				if ((*elem.n6).updated == false){								// If the node hasn't been updated yet...
					(*elem.n6).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n6).collfib.y = coll_fib.y;
					(*elem.n6).collfib.z = coll_fib.z;
					(*elem.n6).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n6).updated = true;}										// Set the updated flag
				else if ((*elem.n6).updated == true){							// If the node has been updated...
					(*elem.n6).collfib.x = ((*elem.n6).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n6).collfib.y = ((*elem.n6).collfib.y + coll_fib.y)/2;
					(*elem.n6).collfib.z = ((*elem.n6).collfib.z + coll_fib.z)/2;
					(*elem.n6).ecm_den = ((*elem.n6).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 7){													// Upper back left node
				if ((*elem.n7).updated == false){								// If the node hasn't been updated yet...
					(*elem.n7).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n7).collfib.y = coll_fib.y;
					(*elem.n7).collfib.z = coll_fib.z;
					(*elem.n7).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n7).updated = true;}										// Set the updated flag
				else if ((*elem.n7).updated == true){							// If the node has been updated...
					(*elem.n7).collfib.x = ((*elem.n7).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n7).collfib.y = ((*elem.n7).collfib.y + coll_fib.y)/2;
					(*elem.n7).collfib.z = ((*elem.n7).collfib.z + coll_fib.z)/2;
					(*elem.n7).ecm_den = ((*elem.n7).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 8){													// Upper back right node
				if ((*elem.n8).updated == false){								// If the node hasn't been updated yet...
					(*elem.n8).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n8).collfib.y = coll_fib.y;
					(*elem.n8).collfib.z = coll_fib.z;
					(*elem.n8).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n8).updated = true;}										// Set the updated flag
				else if ((*elem.n8).updated == true){							// If the node has been updated...
					(*elem.n8).collfib.x = ((*elem.n8).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n8).collfib.y = ((*elem.n8).collfib.y + coll_fib.y)/2;
					(*elem.n8).collfib.z = ((*elem.n8).collfib.z + coll_fib.z)/2;
					(*elem.n8).ecm_den = ((*elem.n8).ecm_den + ecm_den)/2;}}		// Average together the new matrix density

		}
	}

	
	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i){								// Turn off the updated flag for all nodes
		grid.nodes[i].updated = false;
		grid.nodes[i].ecm_den_store.push_back(grid.nodes[i].ecm_den);
		grid.nodes[i].ecm_fibril_store.push_back(grid.nodes[i].collfib);}
	
	//update_ecm_den_grad();										  // Update the ECM density gradient based on the solution from FEBio
		
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - calculate_deform_tensor
//		Calculate the deformation gradient tensor
///////////////////////////////////////////////////////////////////////

mat3d FEAngio::calculate_deform_tensor(Elem elem, const double ex, const double ey, const double ez)
{
	mat3d F;															// Deformation gradient tensor
	mat3d dXde;														// The tensor dX/de (derivative of reference position with respect to natural coordinates) 
	
	// Calculate the derviative of the shape functions evaluate at each node
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;	
   	dN1 = shapefun_d1(ex, ey, ez, 1);								
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);

	// Calculate dX/de
	dXde[0][0] = ((*elem.n1).r0.x)*dN1.x + ((*elem.n2).r0.x)*dN2.x + ((*elem.n3).r0.x)*dN3.x + ((*elem.n4).r0.x)*dN4.x + ((*elem.n5).r0.x)*dN5.x + ((*elem.n6).r0.x)*dN6.x + ((*elem.n7).r0.x)*dN7.x + ((*elem.n8).r0.x)*dN8.x;
	dXde[0][1] = ((*elem.n1).r0.x)*dN1.y + ((*elem.n2).r0.x)*dN2.y + ((*elem.n3).r0.x)*dN3.y + ((*elem.n4).r0.x)*dN4.y + ((*elem.n5).r0.x)*dN5.y + ((*elem.n6).r0.x)*dN6.y + ((*elem.n7).r0.x)*dN7.y + ((*elem.n8).r0.x)*dN8.y;
	dXde[0][2] = ((*elem.n1).r0.x)*dN1.z + ((*elem.n2).r0.x)*dN2.z + ((*elem.n3).r0.x)*dN3.z + ((*elem.n4).r0.x)*dN4.z + ((*elem.n5).r0.x)*dN5.z + ((*elem.n6).r0.x)*dN6.z + ((*elem.n7).r0.x)*dN7.z + ((*elem.n8).r0.x)*dN8.z;

	dXde[1][0] = ((*elem.n1).r0.y)*dN1.x + ((*elem.n2).r0.y)*dN2.x + ((*elem.n3).r0.y)*dN3.x + ((*elem.n4).r0.y)*dN4.x + ((*elem.n5).r0.y)*dN5.x + ((*elem.n6).r0.y)*dN6.x + ((*elem.n7).r0.y)*dN7.x + ((*elem.n8).r0.y)*dN8.x;
	dXde[1][1] = ((*elem.n1).r0.y)*dN1.y + ((*elem.n2).r0.y)*dN2.y + ((*elem.n3).r0.y)*dN3.y + ((*elem.n4).r0.y)*dN4.y + ((*elem.n5).r0.y)*dN5.y + ((*elem.n6).r0.y)*dN6.y + ((*elem.n7).r0.y)*dN7.y + ((*elem.n8).r0.y)*dN8.y;
	dXde[1][2] = ((*elem.n1).r0.y)*dN1.z + ((*elem.n2).r0.y)*dN2.z + ((*elem.n3).r0.y)*dN3.z + ((*elem.n4).r0.y)*dN4.z + ((*elem.n5).r0.y)*dN5.z + ((*elem.n6).r0.y)*dN6.z + ((*elem.n7).r0.y)*dN7.z + ((*elem.n8).r0.y)*dN8.z;

	dXde[2][0] = ((*elem.n1).r0.z)*dN1.x + ((*elem.n2).r0.z)*dN2.x + ((*elem.n3).r0.z)*dN3.x + ((*elem.n4).r0.z)*dN4.x + ((*elem.n5).r0.z)*dN5.x + ((*elem.n6).r0.z)*dN6.x + ((*elem.n7).r0.z)*dN7.x + ((*elem.n8).r0.z)*dN8.x;
	dXde[2][1] = ((*elem.n1).r0.z)*dN1.y + ((*elem.n2).r0.z)*dN2.y + ((*elem.n3).r0.z)*dN3.y + ((*elem.n4).r0.z)*dN4.y + ((*elem.n5).r0.z)*dN5.y + ((*elem.n6).r0.z)*dN6.y + ((*elem.n7).r0.z)*dN7.y + ((*elem.n8).r0.z)*dN8.y;
	dXde[2][2] = ((*elem.n1).r0.z)*dN1.z + ((*elem.n2).r0.z)*dN2.z + ((*elem.n3).r0.z)*dN3.z + ((*elem.n4).r0.z)*dN4.z + ((*elem.n5).r0.z)*dN5.z + ((*elem.n6).r0.z)*dN6.z + ((*elem.n7).r0.z)*dN7.z + ((*elem.n8).r0.z)*dN8.z;

	// Calculate the tensor dM, which is (dX/de)^-T * dN
	vec3d dM1; vec3d dM2; vec3d dM3; vec3d dM4; vec3d dM5; vec3d dM6; vec3d dM7; vec3d dM8;
	
	mat3d dXde_inv_trans;
	dXde_inv_trans = (dXde.inverse()).transpose();

	dM1 = dXde_inv_trans*dN1;
	dM2 = dXde_inv_trans*dN2;
	dM3 = dXde_inv_trans*dN3;
	dM4 = dXde_inv_trans*dN4;
	dM5 = dXde_inv_trans*dN5;
	dM6 = dXde_inv_trans*dN6;
	dM7 = dXde_inv_trans*dN7;
	dM8 = dXde_inv_trans*dN8;

	// Calculate F 
	F[0][0] = ((*elem.n1).rt.x)*dM1.x + ((*elem.n2).rt.x)*dM2.x + ((*elem.n3).rt.x)*dM3.x + ((*elem.n4).rt.x)*dM4.x + ((*elem.n5).rt.x)*dM5.x + ((*elem.n6).rt.x)*dM6.x + ((*elem.n7).rt.x)*dM7.x + ((*elem.n8).rt.x)*dM8.x;
	F[0][1] = ((*elem.n1).rt.x)*dM1.y + ((*elem.n2).rt.x)*dM2.y + ((*elem.n3).rt.x)*dM3.y + ((*elem.n4).rt.x)*dM4.y + ((*elem.n5).rt.x)*dM5.y + ((*elem.n6).rt.x)*dM6.y + ((*elem.n7).rt.x)*dM7.y + ((*elem.n8).rt.x)*dM8.y;
	F[0][2] = ((*elem.n1).rt.x)*dM1.z + ((*elem.n2).rt.x)*dM2.z + ((*elem.n3).rt.x)*dM3.z + ((*elem.n4).rt.x)*dM4.z + ((*elem.n5).rt.x)*dM5.z + ((*elem.n6).rt.x)*dM6.z + ((*elem.n7).rt.x)*dM7.z + ((*elem.n8).rt.x)*dM8.z;

	F[1][0] = ((*elem.n1).rt.y)*dM1.x + ((*elem.n2).rt.y)*dM2.x + ((*elem.n3).rt.y)*dM3.x + ((*elem.n4).rt.y)*dM4.x + ((*elem.n5).rt.y)*dM5.x + ((*elem.n6).rt.y)*dM6.x + ((*elem.n7).rt.y)*dM7.x + ((*elem.n8).rt.y)*dM8.x;
	F[1][1] = ((*elem.n1).rt.y)*dM1.y + ((*elem.n2).rt.y)*dM2.y + ((*elem.n3).rt.y)*dM3.y + ((*elem.n4).rt.y)*dM4.y + ((*elem.n5).rt.y)*dM5.y + ((*elem.n6).rt.y)*dM6.y + ((*elem.n7).rt.y)*dM7.y + ((*elem.n8).rt.y)*dM8.y;
	F[1][2] = ((*elem.n1).rt.y)*dM1.z + ((*elem.n2).rt.y)*dM2.z + ((*elem.n3).rt.y)*dM3.z + ((*elem.n4).rt.y)*dM4.z + ((*elem.n5).rt.y)*dM5.z + ((*elem.n6).rt.y)*dM6.z + ((*elem.n7).rt.y)*dM7.z + ((*elem.n8).rt.y)*dM8.z;

	F[2][0] = ((*elem.n1).rt.z)*dM1.x + ((*elem.n2).rt.z)*dM2.x + ((*elem.n3).rt.z)*dM3.x + ((*elem.n4).rt.z)*dM4.x + ((*elem.n5).rt.z)*dM5.x + ((*elem.n6).rt.z)*dM6.x + ((*elem.n7).rt.z)*dM7.x + ((*elem.n8).rt.z)*dM8.x;
	F[2][1] = ((*elem.n1).rt.z)*dM1.y + ((*elem.n2).rt.z)*dM2.y + ((*elem.n3).rt.z)*dM3.y + ((*elem.n4).rt.z)*dM4.y + ((*elem.n5).rt.z)*dM5.y + ((*elem.n6).rt.z)*dM6.y + ((*elem.n7).rt.z)*dM7.y + ((*elem.n8).rt.z)*dM8.y;
	F[2][2] = ((*elem.n1).rt.z)*dM1.z + ((*elem.n2).rt.z)*dM2.z + ((*elem.n3).rt.z)*dM3.z + ((*elem.n4).rt.z)*dM4.z + ((*elem.n5).rt.z)*dM5.z + ((*elem.n6).rt.z)*dM6.z + ((*elem.n7).rt.z)*dM7.z + ((*elem.n8).rt.z)*dM8.z;

	return F;
}



///////////////////////////////////////////////////////////////////////
// FEAngio - shapefun_d1
//		Evaluate the gradient of the shape functions
///////////////////////////////////////////////////////////////////////

vec3d FEAngio::shapefun_d1(const double xix, const double xiy, const double xiz, int node)
{
    vec3d out;														// Output vector
	
	if (node == 1)													// dN1/de
    {
        out.x = -(1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 - xix)*(1 - xiz)/8.;
		out.z = -(1 - xix)*(1 - xiy)/8.;
    }
    
    if (node == 2)													// dN2/de
    {
        out.x = (1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 - xiy)/8.;
    }    
    
    if (node == 3)													// dN3/de
    {
        out.x = -(1 + xiy)*(1 - xiz)/8.;
        out.y = (1 - xix)*(1 - xiz)/8.;
        out.z = -(1 - xix)*(1 + xiy)/8.;
    }    
    
    if (node == 4)													// dN4/de
    {
        out.x = (1 + xiy)*(1 - xiz)/8.;
        out.y = (1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 + xiy)/8.;
    }    
    
    if (node == 5)													// dN5/de
    {
        out.x = -(1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 - xiy)/8.;
    }    

    if (node == 6)													// dN6/de
    {
        out.x = (1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 - xiy)/8.;
    }    
        
    if (node == 7)													// dN7/de
    {
        out.x = -(1 + xiy)*(1 + xiz)/8.;
        out.y = (1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 + xiy)/8.;
    }   

    if (node == 8)													// dN8/de
    {
        out.x = (1 + xiy)*(1 + xiz)/8.;
        out.y = (1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 + xiy)/8.;
    }   

    return out;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////

void FEAngio::adjust_mesh_stiffness(FEModel& fem)
{
	if (comp_mat == 0)													// If a composite consitutive model isn't being used, exit
		return;
		
	int elem_num = 0;													// Element number
	list<Segment>::iterator frag_it;									// Iterator for the segment container FRAG
	vec3d vess_vect;													// Vessel vector

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i)
	{									
		grid.ebin[i].alpha = 0.;											// Set the vessel volume fraction, alpha, to zero
		grid.ebin[i].fiber_orient.x = 0.;									// Set the vessel orientation vector to 0 (this is the element fiber_orient vector, which does not contain collagen fiber orientation information but rather the material fiber direction for a transversely isotropic material model)
		grid.ebin[i].fiber_orient.y = 0.;
		grid.ebin[i].fiber_orient.z = 0.;
	}
	
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1/(double)Nsub;									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

	for (frag_it = cult.m_frag.begin(); frag_it != cult.m_frag.end(); ++frag_it)		// For each segment...
	{
		Segment seg;													// Segment placeholder
		Segment subunit;												// Segment subdivision placeholder
	
		seg = (*frag_it);												// Obtain the segment
		
		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length > 0.){											// If it's a +1 segment...
					subunit.m_tip[0].rt = seg.m_tip[0].rt;
				}
				else{															// If it's a -1 segment...
					subunit.m_tip[1].rt = seg.m_tip[1].rt;
					subunit.length = -1.;
				}
		}
			
			// Calculate the subdivision
			if (seg.length > 0.){										// If it's a +1 segment...			
				subunit.m_tip[1].rt = subunit.m_tip[0].rt + seg.uvect*(sub_scale*seg.length);     
			}
			else{														// If it's a -1 segment...
				subunit.m_tip[0].rt = subunit.m_tip[1].rt + seg.uvect*(sub_scale*seg.length);     
			}

			subunit.findlength();										// Find the length of the subdivision
			
			mid = (subunit.m_tip[1].rt + subunit.m_tip[0].rt)*0.5;

			elem_num = grid.findelem(mid.x, mid.y, mid.z);				// Find the element that the midpoint is within

			// Calculate the orientation of the subdivision
			if (seg.length > 0.){										// If it's a +1 segment...
				vess_vect = subunit.m_tip[1].rt - subunit.m_tip[0].rt;
			}
			else{														// If it's a -1 segment...
				vess_vect = subunit.m_tip[0].rt - subunit.m_tip[1].rt;
			}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
				vess_vect = vess_vect/vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				elem_volume = grid.ebin[elem_num].volume;					// Calculate the volume of the element
				subunit_volume = pi*(m_vessel_width/2.)*(m_vessel_width/2.)*fabs(subunit.length);		// Find the volume of the subdivision
				volume_fraction = subunit_volume/elem_volume;				// Calculate the volume fraction

				grid.ebin[elem_num].alpha = grid.ebin[elem_num].alpha + volume_fraction;	// Add the volume fraction for each subdivision to alpha
			
				// Calculate the vessel orientation vector 
				if ((grid.ebin[elem_num].fiber_orient.x == 0) && (grid.ebin[elem_num].fiber_orient.y == 0) && (grid.ebin[elem_num].fiber_orient.z == 0)){	// If the vessel orientation vector hasn't been assigned yet...
					grid.ebin[elem_num].fiber_orient.x = vess_vect.x;			// Set the vessel orientation vector					
					grid.ebin[elem_num].fiber_orient.y = vess_vect.y;
					grid.ebin[elem_num].fiber_orient.z = vess_vect.z;}
				else{														// If it has been...	
					grid.ebin[elem_num].fiber_orient.x = (grid.ebin[elem_num].fiber_orient.x + vess_vect.x)/2;	// Average together the vessel orientation vector
					grid.ebin[elem_num].fiber_orient.y = (grid.ebin[elem_num].fiber_orient.y + vess_vect.y)/2;
					grid.ebin[elem_num].fiber_orient.z = (grid.ebin[elem_num].fiber_orient.z + vess_vect.z)/2;}
			}
			
			// Set the origin of the next subdivision to the end of the current one
			if (seg.length > 0.){
				subunit.m_tip[0].rt = subunit.m_tip[1].rt;
			}
			else{
				subunit.m_tip[1].rt = subunit.m_tip[0].rt;
			}
		}
	}
	
	double alpha = 0.;													// Volume fraction for the composite material model
	vec3d e1; vec3d e2; vec3d e3;										// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
	
	FEMesh& mesh = fem.GetMesh();										// Get the FE mesh
	int J = mesh.Domains();												// Find the number of domains within the mesh
    int num_elem = 0;

	for (int k = 0; k < J; ++k){
		FEDomain& d = mesh.Domain(k);										// Obtain the domain
		
		for (int j = 0; j < d.Elements(); ++j)								// For each element within the domain...
		{
			FEElement& e = d.ElementRef(j);										// Obtain the element from the domain
			int nint = e.GaussPoints();											// Obtain the number of gauss points
		
			alpha = grid.ebin[num_elem].alpha;											// Obtain alpha from the grid element

			e1.x = grid.ebin[num_elem].fiber_orient.x;									// Set e1 to the vessel orientation vector
			e1.y = grid.ebin[num_elem].fiber_orient.y;
			e1.z = grid.ebin[num_elem].fiber_orient.z;

			if ((e1.x == 0) && (e1.y == 0) && (e1.z == 0)){						// If there is not vessels in the element, set the material basis to the global coordinate basis
				e1 = vec3d(1,0,0);
				e2 = vec3d(0,1,0);
				e3 = vec3d(0,0,1);}
			else{																// Else, set the other two directions to be orthogonal to the vessel orientation
				e2.y = 1;
				e2 = e1^e2;
				e3 = e1^e2;}
		
			for (int n = 0; n < nint; ++n)										// For each gauss point...
			{
				//FEMaterialPoint& mp = *e.m_State[n];								// Obtain the material point
				//FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();	// Obtain the mixture material point
				//pt.m_w[0] = alpha;													// Set the first weight factor to alpha
				//pt.m_w[1] = 1.0 - alpha;											// Set the second weight factor to 1 - alpha
			
				FEMaterialPoint& mp = *(e.GetMaterialPoint(n)->GetPointData(1));    // returns the second component of the mixture
				FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>(); // get the mixture material point
				pt.m_w[0] = alpha;
				pt.m_w[1] = 1.0 - alpha;

				if (comp_mat == 2){													// If the transversely isotropic material is being used...
					FEElasticMaterialPoint& pt2 = *mp.ExtractData<FEElasticMaterialPoint>();
					pt2.m_Q[0][0] = e1.x;												// Set the first column of Q to e1
					pt2.m_Q[1][0] = e1.y;
					pt2.m_Q[2][0] = e1.z;
					pt2.m_Q[0][1] = e2.x;												// Set the second column of Q to e2
					pt2.m_Q[1][1] = e2.y;
					pt2.m_Q[2][1] = e2.z;
					pt2.m_Q[0][2] = e3.x;												// Set the third column of Q to e3
					pt2.m_Q[1][2] = e3.y;
					pt2.m_Q[2][2] = e3.z;}
				}

			num_elem++;
		}
	}

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - initialize_grid_volume
//		Set the initial grid volume
///////////////////////////////////////////////////////////////////////

void FEAngio::initialize_grid_volume()
{
	mat3d Jacob_mat;														// Jacobian matrix
	double ex = 0.; double ey = 0.; double ez = 0.;						// Position of the centroid in natural coordinates
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;	// Arrays containing the derivatives of the shape functions
	double volume = 0.;													// Element volume

	// Calculate the value of the derviatives of the shape functions evaulated at the centroid
	dN1 = shapefun_d1(ex, ey, ez, 1);
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);
		
	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){									// For each element within the mesh...
		Elem& elem = grid.ebin[i];												// Obtain the element
		
		// Construct the Jacobian matrix
		Jacob_mat[0][0] = ((*elem.n1).rt.x)*dN1.x + ((*elem.n2).rt.x)*dN2.x + ((*elem.n3).rt.x)*dN3.x + ((*elem.n4).rt.x)*dN4.x + ((*elem.n5).rt.x)*dN5.x + ((*elem.n6).rt.x)*dN6.x + ((*elem.n7).rt.x)*dN7.x + ((*elem.n8).rt.x)*dN8.x;
		Jacob_mat[0][1] = ((*elem.n1).rt.x)*dN1.y + ((*elem.n2).rt.x)*dN2.y + ((*elem.n3).rt.x)*dN3.y + ((*elem.n4).rt.x)*dN4.y + ((*elem.n5).rt.x)*dN5.y + ((*elem.n6).rt.x)*dN6.y + ((*elem.n7).rt.x)*dN7.y + ((*elem.n8).rt.x)*dN8.y;
		Jacob_mat[0][2] = ((*elem.n1).rt.x)*dN1.z + ((*elem.n2).rt.x)*dN2.z + ((*elem.n3).rt.x)*dN3.z + ((*elem.n4).rt.x)*dN4.z + ((*elem.n5).rt.x)*dN5.z + ((*elem.n6).rt.x)*dN6.z + ((*elem.n7).rt.x)*dN7.z + ((*elem.n8).rt.x)*dN8.z;
		
		Jacob_mat[1][0] = ((*elem.n1).rt.y)*dN1.x + ((*elem.n2).rt.y)*dN2.x + ((*elem.n3).rt.y)*dN3.x + ((*elem.n4).rt.y)*dN4.x + ((*elem.n5).rt.y)*dN5.x + ((*elem.n6).rt.y)*dN6.x + ((*elem.n7).rt.y)*dN7.x + ((*elem.n8).rt.y)*dN8.x;
		Jacob_mat[1][1] = ((*elem.n1).rt.y)*dN1.y + ((*elem.n2).rt.y)*dN2.y + ((*elem.n3).rt.y)*dN3.y + ((*elem.n4).rt.y)*dN4.y + ((*elem.n5).rt.y)*dN5.y + ((*elem.n6).rt.y)*dN6.y + ((*elem.n7).rt.y)*dN7.y + ((*elem.n8).rt.y)*dN8.y;
		Jacob_mat[1][2] = ((*elem.n1).rt.y)*dN1.z + ((*elem.n2).rt.y)*dN2.z + ((*elem.n3).rt.y)*dN3.z + ((*elem.n4).rt.y)*dN4.z + ((*elem.n5).rt.y)*dN5.z + ((*elem.n6).rt.y)*dN6.z + ((*elem.n7).rt.y)*dN7.z + ((*elem.n8).rt.y)*dN8.z;
		
		Jacob_mat[2][0] = ((*elem.n1).rt.z)*dN1.x + ((*elem.n2).rt.z)*dN2.x + ((*elem.n3).rt.z)*dN3.x + ((*elem.n4).rt.z)*dN4.x + ((*elem.n5).rt.z)*dN5.x + ((*elem.n6).rt.z)*dN6.x + ((*elem.n7).rt.z)*dN7.x + ((*elem.n8).rt.z)*dN8.x;
		Jacob_mat[2][1] = ((*elem.n1).rt.z)*dN1.y + ((*elem.n2).rt.z)*dN2.y + ((*elem.n3).rt.z)*dN3.y + ((*elem.n4).rt.z)*dN4.y + ((*elem.n5).rt.z)*dN5.y + ((*elem.n6).rt.z)*dN6.y + ((*elem.n7).rt.z)*dN7.y + ((*elem.n8).rt.z)*dN8.y;
		Jacob_mat[2][2] = ((*elem.n1).rt.z)*dN1.z + ((*elem.n2).rt.z)*dN2.z + ((*elem.n3).rt.z)*dN3.z + ((*elem.n4).rt.z)*dN4.z + ((*elem.n5).rt.z)*dN5.z + ((*elem.n6).rt.z)*dN6.z + ((*elem.n7).rt.z)*dN7.z + ((*elem.n8).rt.z)*dN8.z;

		if (Jacob_mat.det() != 0.)											// Calculate the volume from the Jacobian matrix using the determinant
			volume = 8*Jacob_mat.det();

		grid.ebin[i].volume = volume;										// Set the element volume
		grid.ebin[i].volume0 = volume;}										// Set the element initial volume
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_grid_volume
//		Update grid volume for each element after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_grid_volume()
{
	mat3d F;																// Deformation gradient tensor
	double Jacob = 0.;													// Jacobian (i.e., determinant of F)
	double ex = 0.; double ey = 0.; double ez = 0.;						// Position of the centroid in natural coordinates
	double new_volume = 0.;												// New element volume

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){									// For each element within the mesh...
		Elem& elem = grid.ebin[i];												// Obtain the element
		F = calculate_deform_tensor(elem, ex, ey, ez);						// Calculate the deformation gradient tensor
		Jacob = F.det();													// Calculate the Jacobian by taking the determinant of F

		new_volume = Jacob*elem.volume0;									// Calculate the new element volume using the Jacobian
		grid.ebin[i].volume = new_volume;}									// Store the new element volume
	
	return;
}

///////////////////////////////////////////////////////////////////////
// FEAngio - output_params
//		Output parameter values after the simulation ends
///////////////////////////////////////////////////////////////////////

void FEAngio::output_params()
{
	FILE *param_stream;													// Parameter output file stream                                                                                                                     
	param_stream = fopen("out_params.ang","wt");						// Output the parameter output file
	
	fprintf(param_stream,"a = %5.5f \n",m_sproutf);						// Print the sprout force magnitude
	fprintf(param_stream,"tip range = %5.5f \n",m_tip_range);				// Print the sprout force range
	fprintf(param_stream,"phi_stiff_factor = %5.5f \n",phi_stiff_factor);	// Print the displacement stiffness factor
	fprintf(param_stream,"total_body_force = %10.5i \n",total_bdyf);		// Print the total number of body forces
	fclose(param_stream);												// Close the parameter output file

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - save_cropped_vessels
//		Output parameter values after the simulation ends
///////////////////////////////////////////////////////////////////////

void FEAngio::save_cropped_vessels()
{
	FILE *data_stream;                                                          // Open stream to 'grid.ang' (stream2)                                                                   
	data_stream = fopen("out_cropped_data.ang","wt");
	
	FILE *seg_conn_stream;                                                           // Open stream to 'out_seg_conn.ang' (stream)
	seg_conn_stream = fopen("out_cropped_seg_conn.ang","wt");

	list<Segment>::iterator it;
		
	double xmin = 0.; double xmax = 0.; double ymin = 0.; double ymax = 0.; double zmin = 0.; double zmax = 0.;
	
	xmin = ((grid.xrange[1] + grid.xrange[0])/2) - 2548/2;
	xmax = ((grid.xrange[1] + grid.xrange[0])/2) + 2548/2;
	ymin = ((grid.yrange[1] + grid.yrange[0])/2) - 2548/2;
	ymax = ((grid.yrange[1] + grid.yrange[0])/2) + 2548/2;
	zmin = ((grid.zrange[1] + grid.zrange[0])/2) - 1500/2;
	zmax = ((grid.zrange[1] + grid.zrange[0])/2) + 1500/2;

	for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		vec3d& r0 = it->m_tip[0].rt;
		vec3d& r1 = it->m_tip[1].rt;
		if (((r0.x <= xmax) && (r0.x >= xmin)) && ((r1.x <= xmax) && (r1.x >= xmin))){
			if (((r0.y <= ymax) && (r0.y >= ymin)) && ((r1.y <= ymax) && (r1.y >= ymin))){
				if (((r0.z <= zmax) && (r0.z >= zmin)) && ((r1.z <= zmax) && (r1.z >= zmin))){
					fprintf(data_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->seg_num,it->label);  // Write to data.ang		
					fprintf(seg_conn_stream,"%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i\n",it->seg_num,it->seg_conn[0][0],it->seg_conn[0][1],it->seg_conn[1][0],it->seg_conn[1][1]);  // Write to seg_conn.ang
				}
			}
		}
	}
	
	fclose(data_stream);
	fclose(seg_conn_stream);

	return;

}

///////////////////////////////////////////////////////////////////////
// FEAngio - update_ecm_den_grad
//		Calculate the density gradient for each element (this function may not work properly)
///////////////////////////////////////////////////////////////////////

void FEAngio::update_ecm_den_grad()
{
	vec3d elem_den_grad;
	double ex = 0.; double ey = 0.; double ez = 0.;	
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;
	
	dN1 = shapefun_d1(ex, ey, ez, 1);
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);
	
	int NN = grid.Nodes();
	for (int j = 0; j < NN; ++j){
		grid.nodes[j].updated = false;}


	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){
		Elem& elem = grid.ebin[i];
		
		elem_den_grad.x = ((*elem.n1).ecm_den)*dN1.x + ((*elem.n2).ecm_den)*dN2.x + ((*elem.n3).ecm_den)*dN3.x + ((*elem.n4).ecm_den)*dN4.x + ((*elem.n5).ecm_den)*dN5.x + ((*elem.n6).ecm_den)*dN6.x + ((*elem.n7).ecm_den)*dN7.x + ((*elem.n8).ecm_den)*dN8.x;
		elem_den_grad.y = ((*elem.n1).ecm_den)*dN1.y + ((*elem.n2).ecm_den)*dN2.y + ((*elem.n3).ecm_den)*dN3.y + ((*elem.n4).ecm_den)*dN4.y + ((*elem.n5).ecm_den)*dN5.y + ((*elem.n6).ecm_den)*dN6.y + ((*elem.n7).ecm_den)*dN7.y + ((*elem.n8).ecm_den)*dN8.y;
		elem_den_grad.z = ((*elem.n1).ecm_den)*dN1.z + ((*elem.n2).ecm_den)*dN2.z + ((*elem.n3).ecm_den)*dN3.z + ((*elem.n4).ecm_den)*dN4.z + ((*elem.n5).ecm_den)*dN5.z + ((*elem.n6).ecm_den)*dN6.z + ((*elem.n7).ecm_den)*dN7.z + ((*elem.n8).ecm_den)*dN8.z;
	
		if ((*elem.n1).updated == false){
			(*elem.n1).ecm_den_grad = elem_den_grad;
			(*elem.n1).updated = true;}
		else if ((*elem.n1).updated == true){
			(*elem.n1).ecm_den_grad = (elem_den_grad + (*elem.n1).ecm_den_grad)*0.5;}

		if ((*elem.n2).updated == false){
			(*elem.n2).ecm_den_grad = elem_den_grad;
			(*elem.n2).updated = true;}
		else if ((*elem.n2).updated == true){
			(*elem.n2).ecm_den_grad = (elem_den_grad + (*elem.n2).ecm_den_grad)*0.5;}
		
		if ((*elem.n3).updated == false){
			(*elem.n3).ecm_den_grad = elem_den_grad;
			(*elem.n3).updated = true;}
		else if ((*elem.n3).updated == true){
			(*elem.n3).ecm_den_grad = (elem_den_grad + (*elem.n3).ecm_den_grad)*0.5;}

		if ((*elem.n4).updated == false){
			(*elem.n4).ecm_den_grad = elem_den_grad;
			(*elem.n4).updated = true;}
		else if ((*elem.n4).updated == true){
			(*elem.n4).ecm_den_grad = (elem_den_grad + (*elem.n4).ecm_den_grad)*0.5;}

		if ((*elem.n5).updated == false){
			(*elem.n5).ecm_den_grad = elem_den_grad;
			(*elem.n5).updated = true;}
		else if ((*elem.n5).updated == true){
			(*elem.n5).ecm_den_grad = (elem_den_grad + (*elem.n5).ecm_den_grad)*0.5;}

		if ((*elem.n6).updated == false){
			(*elem.n6).ecm_den_grad = elem_den_grad;
			(*elem.n6).updated = true;}
		else if ((*elem.n6).updated == true){
			(*elem.n6).ecm_den_grad = (elem_den_grad + (*elem.n6).ecm_den_grad)*0.5;}

		if ((*elem.n7).updated == false){
			(*elem.n7).ecm_den_grad = elem_den_grad;
			(*elem.n7).updated = true;}
		else if ((*elem.n7).updated == true){
			(*elem.n7).ecm_den_grad = (elem_den_grad + (*elem.n7).ecm_den_grad)*0.5;}

		if ((*elem.n8).updated == false){
			(*elem.n8).ecm_den_grad = elem_den_grad;
			(*elem.n8).updated = true;}
		else if ((*elem.n8).updated == true){
			(*elem.n8).ecm_den_grad = (elem_den_grad + (*elem.n8).ecm_den_grad)*0.5;}
	}
	
	return;
}


void FEAngio::update_sprout_stress_scaling()
{
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	if (m_pmat)
		m_pmat->scale = y0 + a/(1 + exp(-(m_t - x0)/b));
	
	return;
}

void FEAngio::circ_gel()
{
	double xmax = grid.xrange[1];
	double ymax = grid.yrange[1];

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){
		Elem& elem = grid.ebin[i];

		// Check face 1 for right boundary
		if ((((*elem.n1).rt.x == xmax) && ((*elem.n2).rt.x == xmax) && ((*elem.n5).rt.x == xmax) && ((*elem.n6).rt.x == xmax)) || (((*elem.n1).rt.y == ymax) && ((*elem.n2).rt.y == ymax) && ((*elem.n5).rt.y == ymax) && ((*elem.n6).rt.y == ymax))){
			grid.ebin[i].f1.BC = true;
			grid.ebin[i].f1.bc_type = 'w';
		}

		// Check face 2 for right boundary
		if ((((*elem.n2).rt.x == xmax) && ((*elem.n4).rt.x == xmax) && ((*elem.n6).rt.x == xmax) && ((*elem.n8).rt.x == xmax)) || (((*elem.n2).rt.y == ymax) && ((*elem.n4).rt.y == ymax) && ((*elem.n6).rt.y == ymax) && ((*elem.n8).rt.y == ymax))){
			grid.ebin[i].f2.BC = true;
			grid.ebin[i].f2.bc_type = 'w';
		}

		// Check face 3 for right boundary
		if ((((*elem.n3).rt.x == xmax) && ((*elem.n4).rt.x == xmax) && ((*elem.n7).rt.x == xmax) && ((*elem.n8).rt.x == xmax)) || (((*elem.n3).rt.y == ymax) && ((*elem.n4).rt.y == ymax) && ((*elem.n7).rt.y == ymax) && ((*elem.n8).rt.y == ymax))){
			grid.ebin[i].f3.BC = true;
			grid.ebin[i].f3.bc_type = 'w';
		}

		// Check face 4 for right boundary
		if ((((*elem.n1).rt.x == xmax) && ((*elem.n3).rt.x == xmax) && ((*elem.n5).rt.x == xmax) && ((*elem.n7).rt.x == xmax)) || (((*elem.n1).rt.y == ymax) && ((*elem.n3).rt.y == ymax) && ((*elem.n5).rt.y == ymax) && ((*elem.n7).rt.y == ymax))){
			grid.ebin[i].f4.BC = true;
			grid.ebin[i].f4.bc_type = 'w';
		}

		double xpt; double ypt; 	
		vec3d r;
		double d;
		double ep = 200.;

		//// Check face 1 for circumfrential boundary
		//xpt = ((*elem.n1).x + (*elem.n2).x + (*elem.n5).x + (*elem.n6).x)/4;
		//ypt = ((*elem.n1).y + (*elem.n2).y + (*elem.n5).y + (*elem.n6).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f1.BC = true;
		//	grid.ebin[i].f1.bc_type = 'i';}
	
		// Check face 2 for circumfrential boundary
		xpt = ((*elem.n2).rt.x + (*elem.n4).rt.x + (*elem.n6).rt.x + (*elem.n8).rt.x)/4;
		ypt = ((*elem.n2).rt.y + (*elem.n4).rt.y + (*elem.n6).rt.y + (*elem.n8).rt.y)/4;

		r.x = xpt - xmax; r.y = ypt - ymax;
		d = r.norm();

		if (fabs(d - xmax) <= ep){
			grid.ebin[i].f2.BC = true;
			grid.ebin[i].f2.bc_type = 'i';}

		//// Check face 3 for circumfrential boundary
		//xpt = ((*elem.n3).x + (*elem.n4).x + (*elem.n7).x + (*elem.n8).x)/4;
		//ypt = ((*elem.n3).y + (*elem.n4).y + (*elem.n7).y + (*elem.n8).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f3.BC = true;
		//	grid.ebin[i].f3.bc_type = 'i';}

		//// Check face 4 for circumfrential boundary
		//xpt = ((*elem.n1).x + (*elem.n3).x + (*elem.n5).x + (*elem.n7).x)/4;
		//ypt = ((*elem.n1).y + (*elem.n3).y + (*elem.n5).y + (*elem.n7).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f4.BC = true;
		//	grid.ebin[i].f4.bc_type = 'i';}
	}

	return;
}

///////////////////////////////////////////////////////////////////////
// enfore_fiber_BCS
///////////////////////////////////////////////////////////////////////

void FEAngio::enforce_fiber_BCS(Node &node, bool circ)
{
	if (circ == true){
		vec3d r;
		r.x = node.rt.x - grid.xrange[1];
		r.y = node.rt.y - grid.yrange[1];

		r = r/r.norm();

		mat3d I(1.0);
		mat3d M = r&r;

		vec3d r_new; r_new = (I - M)*node.collfib; r_new = r_new/r_new.norm();

		node.collfib.x = r_new.x; node.collfib.y = r_new.y; node.collfib.z = r_new.z;

		return;
	}
		
	if (m_cgelbc == 'l'){
		//Node on the front face
		if (node.rt.y == grid.yrange[0])
			node.collfib.y = 0;

		// Node on the back face
		//if (node.rt.y == grid.yrange[1])
			//node.collfib.y = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.collfib.z = 0;

		// Node on the top face
		//if (node.rt.z == grid.zrange[1])
			//node.collfib.z = 0;
	}

	if (m_cgelbc == 's'){
		// Node on the right face
		//if (node.rt.x == grid.xrange[1])
			//node.collfib.x = 0;

		// Node on the left face
		if (node.rt.x == grid.xrange[0])
			node.collfib.x = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.collfib.z = 0;

		// Node on the top face
		//if (node.rt.z == grid.zrange[1])
			//node.collfib.z = 0;
	}

	if (m_cgelbc == 'u'){
		//Node on the front face
		if (node.rt.y == grid.yrange[0])
			node.collfib.y = 0;

		// Node on the right face
		//if (node.rt.x == grid.xrange[0])
			//node.collfib.x = 0;

		// Node on the back face
		//if (node.rt.y == grid.yrange[1])
			//node.collfib.y = 0;

		// Node on the left face
		if (node.rt.x == grid.xrange[0])
			node.collfib.x = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.collfib.z = 0;

		// Node on the top face
		//if (node.rt.z == grid.zrange[1])
			//node.collfib.z = 0;
	}

	if (m_cgelbc == 'n'){
	}

	return;
}

//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientation.
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat)
{
	// get the material's coordinate system
	FECoordSysMap* pmap = pmat->GetCoordinateSystemMap();
//	if (pmap == 0) return false;

	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the fiber array
	int N = mesh.Nodes();
	fiber.resize(N, vec3d(0,0,0));

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local fiber orientation at integration points
			vector<double> fx(nint), fy(nint), fz(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMaterialPoint& pt = *mpoint->ExtractData<FEElasticMaterialPoint>();
				mat3d m = pt.m_Q;
				//mat3d m = pmap->LocalElementCoord(el, n);
				
				// grab the first column as the fiber orientation
				fx[n] = m[0][0];
				fy[n] = m[1][0];
				fz[n] = m[2][0];
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln), gy(neln), gz(neln);
			el.project_to_nodes(&fx[0], &gx[0]);
			el.project_to_nodes(&fy[0], &gy[0]);
			el.project_to_nodes(&fz[0], &gz[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				vec3d f(gx[i], gy[i], gz[i]);
				fiber[ni] += f;
			}
		}
	}

	// normalize the fibers
	for (int i=0; i<N; ++i) fiber[i].unit();

	// If we get here, all is well.
	return true;
}

//-----------------------------------------------------------------------------
// create a density map based on material density parameter per point
bool CreateDensityMap(vector<double>& density, FEMaterial* pmat)
{
	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the density array
	int N = mesh.Nodes();
	density.resize(N, 0.0);

	// loop over all domains
	for (int nd = 0; nd < mesh.Domains(); ++nd)
	{
		FESolidDomain& dom = static_cast<FESolidDomain&>(mesh.Domain(nd));
		
		// loop over all elements in the domain
		int NE = dom.Elements();
		for (int ne=0; ne<NE; ++ne)
		{
			FESolidElement& el = dom.Element(ne);
			int neln = el.Nodes();
			int nint = el.GaussPoints();

			// local density at integration points
			vector<double> den(nint);
			for (int n=0; n<nint; ++n)
			{
				// generate a coordinate transformation at this integration point
				FEMaterialPoint* mpoint = el.GetMaterialPoint(n);
				FEElasticMixtureMaterialPoint* mPt = dynamic_cast<FEElasticMixtureMaterialPoint*>(mpoint);

				if(mPt)
				{
					vector<FEMaterialPoint*> mPtV = mPt->m_mp;
					for (int i=0; i<(int)mPtV.size(); ++i)
					{
						FEAngioMaterialPoint* angioPt = dynamic_cast<FEAngioMaterialPoint*>(mPtV[i]);
						if(angioPt)
						{
							den[n] = angioPt->m_D;
							break;
						}
					}
				}
			}

			// now that we have the values at the integration points
			// we need to map it to the nodes
			vector<double> gx(neln);
			el.project_to_nodes(&den[0], &gx[0]);

			// add it to the node accumulators
			for (int i=0; i<neln; ++i)
			{
				int ni = el.m_node[i];
				density[ni] = gx[i];
			}
		}
	}

	// If we get here, all is well.
	return true;
}

///////////////////////////////////////////////////////////////////////
// initBranch
///////////////////////////////////////////////////////////////////////

void FEAngio::initBranch()
{
	list<Segment>::iterator it;
	
	for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)
	{
	     if (float(rand())/RAND_MAX < (m_a1 + m_a2))        // Generate random number between 0 and 1
			it->init_branch = true;                             // If random number less than a1 + a2, then initial branching is true
		else
			it->init_branch = false;                            // If random number is not less than a1 + a2, then initial branching is false
	}
}

//-----------------------------------------------------------------------------
// Update the time step size and current time value.
void FEAngio::updateTime()
{
	m_dt = m_dtA*pow(E,m_dtB/(m_n + m_dtC));
    m_n += 1;

    if (m_dt > (m_maxt-m_t) && (m_maxt - m_t) > 0)       // If dt is bigger than the amount of time left...
		m_dt = m_maxt - m_t;                                       // then just set dt equal to the amount of time left (maxt-t)

	
    m_t = m_t + m_dt;                                    // Update time
}

//-----------------------------------------------------------------------------
// Calculates and stores the total lenght of all vessels.
void FEAngio::updateTotalLength()
{
    m_total_length = 0.;
    list<Segment>::iterator it;
    for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)
    {
        m_total_length += fabs(it->length);    
    }
}


///////////////////////////////////////////////////////////////////////
// Growth
///////////////////////////////////////////////////////////////////////

void FEAngio::Growth()
{
    int k;
    list<Segment>::iterator it;                                 // Declare iterator through list container 'it'
    int test = 0;

	//// Elongate the active vessel tips
    for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)             // Iterate through each vessel fragment (it)
	{
		//it->findphi();
		test = 1;

		for (k=0; k<=1; ++k)                                        // Iterate through each tip of the vessel fragment (k)
		{
	        if (it->m_tip[k].active != 0)                                        // If tip is active (not 0) then create a new segment
			{
		        Segment seg;                                                // Declare SEGMENT object 'seg'
				seg = cult.createNewSeg(it,k);               // CULTURE.createNewSeg(list iterator, grid, data, tip, frag container)
				                                                            // Create new vessel segment on existing segment to 
				                                                            // 'elongate' vessel over time step 
				
				//cult.CheckForIntersection(seg,frag,data,it);                // CULTURE.CheckForIntersection(new segment, fragment container, data, list iterator)
				                                                            // Checks for intersection between new segment and any existing fragment
				                                                               // If true, '3D Insertion' notice will be printed
				
				
				
				cult.m_frag.push_front (seg);                                      // Append new segment at the top of the list 'frag'
				
			}
            
	    }
		    
		Branch(it);												// Determine which segments form new branches

    }

	Fuse();														// Determine which segments form new anastomoses
	
	return;
}



///////////////////////////////////////////////////////////////////////
// updateLength
///////////////////////////////////////////////////////////////////////

void FEAngio::updateLength()
{
	double lc;                                                  // lc - Length calculation obtained from growth function g(t)
    
	lc = m_a/(1.+pow(E,-(m_t-m_x0)/m_b));
	lc -= m_a/(1.+pow(E,-(m_t-m_dt-m_x0)/m_b));
	
	m_vess_length = lc*m_length_adjust;
}


///////////////////////////////////////////////////////////////////////
// Branch
///////////////////////////////////////////////////////////////////////

void FEAngio::Branch(list<Segment>::iterator it)
{
    int k;
    //list<Segment>::iterator it;                                         // Declare iterator through list container 'it'
    
    // Generate a random number between 0 and 1. If that number is less than the branching probability at
    // the current time, or initial branching is activated (data.ini_branch = true), then branching occurs
	
    double den_scale = 1.0;

	vec3d pt = (it->m_tip[1].rt + it->m_tip[0].rt)*0.5;
	
	den_scale = cult.findDenScale(pt.x, pt.y, pt.z);
	
	if ( float(rand())/RAND_MAX < den_scale*m_dt*m_branch_chance/m_t || (it->init_branch == true) )
    {
	    if ((it->BCdead == 0) && (it->anast == 0))                  // Segments that have encountered a boundary condition or formed an anastomoses may not
		{                                                           // form a branch.
	        
			Segment seg;                                            // Declare SEGMENT object 'seg'
			m_num_branches = m_num_branches + 1;              // Iterate the total number of branches +1
			m_branch = true;                                     // Branching flag set to 'true.' This tells the program that the
			                                                        // new vessel segment being created is arising from a branch      
    		it->init_branch = false;
    				
			if (float(rand())/RAND_MAX < 0.5)                       // Randomly determine which node of the parent segment forms the branch
			    k = 0;
		    else
			    k = 1;
    							
			it->m_tip[k].active = sign(0.5f - float(rand())/RAND_MAX);        // Randomly assign the branch to grow as +1 or -1
    				
			seg = cult.createNewSeg(it,k);           // CULTURE.createNewSeg(list iterator, grid, data, tip, frag container)
			                                                            // Create new vessel segment on existing segment to 
			                                                            // create the branching 

            cult.m_num_vessel = cult.m_num_vessel + 1;
		    seg.vessel = cult.m_num_vessel;
    				                    
			it->Recent_branch = 1;
			seg.Recent_branch = 1;                                  // Indicate that this segment just branched, prevents it from 
				                                                        // branching again too soon
    		
			cult.m_frag.push_front (seg);                                  // Append new segment at the top of the list 'frag'
			
			m_branch = false;                                    // Turn off branching flag once branching algorithm is complete
        }
	}                                                               // End Branching
	

	return;
}



///////////////////////////////////////////////////////////////////////
// Fuse
///////////////////////////////////////////////////////////////////////

void FEAngio::Fuse()
{
    list<Segment>::iterator it;
           
    for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		check4anast(it, 0);
        check4anast(it, 1);
	} 
			
	return;
}



///////////////////////////////////////////////////////////////////////
// check4anast
///////////////////////////////////////////////////////////////////////

void FEAngio::check4anast(list<Segment>::iterator it, int k)
{
    if (it->m_tip[k].active == 0)
        return;
	
	if (it->anast == 1)
		return;

	int kk = 0;
    double dist0 = 0.;
    double dist1 = 0.;
    list<Segment>::iterator it2;
	
    for (it2 = cult.m_frag.begin(); it2 != cult.m_frag.end(); ++it2)          // Iterate through all segments in frag list container again (it2)
	{                                                           
	    dist0 = (it->m_tip[k].rt - it2->m_tip[0].rt).norm(); dist0 *= dist0;
		dist1 = (it->m_tip[k].rt - it2->m_tip[1].rt).norm(); dist1 *= dist1;

        anastomose(dist0, dist1, k, it, it2);
    } 
}


///////////////////////////////////////////////////////////////////////
// anastomose
///////////////////////////////////////////////////////////////////////

void FEAngio::anastomose(double dist0, double dist1, int k, list<Segment>::iterator it, list<Segment>::iterator it2)
{
	if ((it->anast == 1) || (it2->anast == 1))
        return;
    
    if (it->label == it2->label)
        return;

    int kk = 9;
    
    if (dist0 <= m_anast_dist)
        kk = 0;
    else if (dist1 <= m_anast_dist)
        kk = 1;
  
    if (kk == 9)
        return;
                                                
	//if (it2->tip[kk] == 0)                                      // Tip-to-tip anastomose only
    //  return;
     
    Segment seg;                                                // Declare SEGMENT object 'seg'
								
	seg = cult.connectSegment(it,it2,k,kk);      // CULTURE.connectSegment(segment 1, segment 2, tip 1, tip 2, grid, data, frag list container)
					                                            // This function will create a segment between to two segments to complete the anastomosis
	cult.m_frag.push_front (seg);                                      // Append new segment at the top of the list 'frag'                  
	it->m_tip[k].active=0;                                               // Deactivate tip of segment 1 after anastomosis
	it2->m_tip[kk].active=0;                                             // Deactivate tip of segment 2 after anastomosis (tip-tip anastomosis only)

    
	return;
}

///////////////////////////////////////////////////////////////////////
// removeErrors
///////////////////////////////////////////////////////////////////////

void FEAngio::removeErrors()
{
	if (kill_off == false){
		double length_limit = m_d;
    
		list<Segment>::iterator it;
    
		for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it){
			if (fabs(it->length) >= 2.0*length_limit){
				it->death_label = -7;
				vec3d& r0 = it->m_tip[0].rt;
				vec3d& r1 = it->m_tip[1].rt;
				fprintf(killed_segs_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->death_label);
				it = cult.m_frag.erase(it);}
			if (it == cult.m_frag.end())
				return;
		}
                
		return;
	}
}

///////////////////////////////////////////////////////////////////////
// find_active_tips()
///////////////////////////////////////////////////////////////////////

void FEAngio::find_active_tips()
{
	list<Segment>::iterator it;
           
	active_tips.clear();
	
	for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		if (it->m_tip[0].active != 0){
			active_tips.push_back(it);}
		else if (it->m_tip[1].active != 0){
			active_tips.push_back(it);}
	} 
			
    return;
}


///////////////////////////////////////////////////////////////////////
// kill_dead_segs
///////////////////////////////////////////////////////////////////////

void FEAngio::kill_dead_segs()
{  
	if (kill_off == false){
		list<Segment>::iterator it;
    
		for (it = cult.m_frag.begin(); it != cult.m_frag.end(); ++it){
			if (it->mark_of_death == true){
				vec3d& r0 = it->m_tip[0].rt;
				vec3d& r1 = it->m_tip[1].rt;
				fprintf(killed_segs_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->death_label);
				it = cult.m_frag.erase(it);}
		}
	}
}
