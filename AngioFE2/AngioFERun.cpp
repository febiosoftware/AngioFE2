#include "StdAfx.h"
#include "AngioFETask.h"
#include "Filein.h"
#include "FEAngio.h"
#include <FECore/FEModel.h>
#include <FECore/log.h>

//-----------------------------------------------------------------------------
// Store a global pointer to the FEAngio object so that other classes
// can access it conveniently (e.g. the plot classes).
// TODO: I'd like to eliminate if possible
FEAngio* pfeangio = 0;

//-----------------------------------------------------------------------------
// Task initialization.
// This allocates the FEAngio object and read the angio input file. It also initializes the 
// FEAngio object as well as the FE model.
bool AngioFETask::Init(const char* inpfile)
{
	// Get the FE model
	FEModel& fem = *GetFEModel();

	// Create the FEAngio class
	pfeangio = new FEAngio(fem);
	FEAngio& feangio = *pfeangio;

	// Read the angio3d input file
	Filein filein;
	if (filein.Input(inpfile, feangio) == false) return false;
	
	// initialize feangio object
	if (feangio.Init() == false) return false;

	// Do the model initialization
	if (fem.Init() == false) return false;

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::FILE_ONLY);

	// all is well
	return true;
}

//-----------------------------------------------------------------------------
bool AngioFETask::Run()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEAngio& feangio = *pfeangio;
	FEAnalysis* pstep = fem.GetCurrentStep();

	feangio.save_vessel_state();									// Output microvessel state
	feangio.save_time();											// Output time information
	feangio.fileout.writeCollFib(feangio.grid, true);				// Output initial collagen fiber orientation

	//// FEBIO - Solve for initial step
	feangio.FE_state++;												// Update the FE state
	fem.Solve();													// Solve the FE model
	feangio.update_grid(mesh);										// Update the position of grid nodes within angio3d using the solution from FEBio
	feangio.update_ECM();											// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)
	feangio.update_grid_volume();									// Update the volume for each element in the grid
	feangio.displace_vessels();										// Displace the microvessels using the solution from FEBio
	pstep->m_dt = pstep->m_dt0;
	feangio.FE_time_step = pstep->m_dt;

	//// ANGIO3D - Apply sprout forces at active growth tips
	feangio.apply_sprout_forces(fem, 1, 0.5);						// Apply sprout forces to the mesh
	feangio.adjust_mesh_stiffness(fem);								// Adjust the stiffness of the mesh based on microvessel volume
	feangio.save_vessel_state();									// Output microvessel state
	//feangio.save_bdy_forces(fem);									// Output body force state
	feangio.save_time();											// Output time

	//// ANGIO3D - Simulate angiogenesis
	feangio.initBranch();											// Handle branching within inital fragments
    feangio.updateTotalLength();									// Update the total vascular length within the simulation   
	
	while (feangio.data.t < feangio.data.maxt)                      // While culture time is less than the maximum time
	{	
		feangio.updateTime();											// Determine the current time step and update time

		feangio.updateLength();											// Determine length of new segments for this growth step        		

		feangio.Growth(fem);											// Growth (Elongation, Branching, and Anastomosis)
    	
		feangio.kill_dead_segs();										// Remove buggy segments
		
		feangio.updateTotalLength();									// Update total vascular length
	    
		feangio.find_active_tips();										// Locate all active growth tips
		
		feangio.adjust_mesh_stiffness(fem);								// Uncomment this if not using a composite material model.  This will update the stiffness of any element that contains microvessel segements 
		
		// Discretize the growth step over N steps in quasi-time to produce smoother results (default N = 2)
		if (feangio.Subgrowth(feangio.m_sub_cycles, fem) == false) break;
		
		feangio.update_ECM();											// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)

		feangio.update_grid_volume();									// Re-calculate the volume of each elment after the deformation using the Jacobian								

		feangio.save_time();											// Output time information	
		
		feangio.fileout.printStatus(feangio);						// Print the status of angio3d to the user    
	}
	
	//// ANGIO3D - Generate output files
	feangio.output_params();										// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	
	//feangio.fileout.printsproutnodes(feangio.sprout_nodes);		// Output the sprout locations
	
	feangio.fileout.dataout(feangio);								// Output data file
		
	feangio.fileout.writeCollFib(feangio.grid, false);				// Output final collagen fiber orientation

	feangio.fileout.writeECMDen(feangio.grid);						// Output final matrix density

	feangio.fileout.writeSegConn(feangio.cult.m_frag);						// Output the segment connectivity data
	
	feangio.fileout.writeECMDenStore(feangio.grid);

	feangio.fileout.writeECMFibrilStore(feangio.grid);

	//feangio.fileout.writeECMDenGrad(feangio.grid);				// Output the ECM density gradient, I don't think this works right...

	feangio.save_cropped_vessels();									// Output the vessels within a selected region of interest

	return true;					
}
