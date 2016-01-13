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
#include <FECore/FEModel.h>
#include "FECore\FESolidDomain.h"
#include "FEBioMech\FEElasticMaterial.h"
#include "FEAngioMaterial.h"
#include "FEBioMech/FEElasticMixture.h"
#include "Elem.h"
#include "angio3d.h"
#include "Culture.h"
#include <iostream>


//-----------------------------------------------------------------------------
// create a map of fiber vectors based on the material's orientatin
bool CreateFiberMap(vector<vec3d>& fiber, FEMaterial* pmat);

// create a density map based on material point density
bool CreateDensityMap(vector<double>& density, FEMaterial* pmat);

//-----------------------------------------------------------------------------
FEAngio::FEAngio(FEModel& fem) : m_fem(fem), m_grid(fem.GetMesh())
{
	// Create the culture
	m_pCult = new Culture(*this);

	// Body force counter
	total_bdyf = 0;
	
	FE_state = 0;

	comp_mat = 0;
	m_pmat = 0;

	phi_stiff_factor = 1.0;
	m_sub_cycles = 2;
    m_ntime = 1;

	// flag for generating fibers (0 = random, 3 = element orientation)
	m_matrix_cond = 0;

	// flatten fiber option (default to false)
	m_bzfibflat = 0;

	// initialize time stepping parameters
	m_time.dt = 0.25;

	// TODO: What are these and make these user parameters
    m_dtA = 0.0637;
    m_dtB = 9.0957;
    m_dtC = 2.6073;

	// vessel_width - Diameter of microvessels (Default: 7 um)
	m_vessel_width = 7;

    m_total_length = 0.;
    
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
}

//-----------------------------------------------------------------------------
FEAngio::~FEAngio()
{
	delete m_pCult;
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

	// Initialize Culture class
	if (m_pCult->Init() == false) return false;

	// create the grid based on the FEBio mesh
	if (m_grid.Init() == false) return false;

	// assign ECM densities to grid nodes
	if (InitECMDensity() == false) return false;

	// assign collagen fibers to grid nodes
	if (InitCollagenFibers() == false) return false;

	// Seed initial fragments 
	// NOTE: must be done after InitECMDensity() and InitCollagenFibers().
	m_pCult->SeedFragments(m_time);

	// Init all the FE stuff
	if (InitFEM() == false) return false;
	
	// start timer
	time(&m_start);

	return true;
}

//-----------------------------------------------------------------------------
// Initialize FE model.
bool FEAngio::InitFEM()
{
	// See if an "angio" material is defined.
	bool bmat = true;
	FEAngioMaterial* pma = FindAngioMaterial(m_fem.GetMaterial(0));
	if (pma == 0) bmat = false;
	else 
	{
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

	//// FEBIO - Initialize the FE model
	m_pmat->SetFEAngio(this);

	// report if the stress or body force approach will be used
	if (bmat)
		felog.printf("Angio material found. Stress approach will be used.");
	else
		felog.printf("Angio materia NOT found. Body-force appraoch will be used.");

	// Do the model initialization
	if (m_fem.Init() == false) return false;

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::FILE_ONLY);

	// --- Output initial state of model ---

	// Output initial microvessel state
	fileout.save_vessel_state(*this);

	// Output time information
	fileout.save_time(*this);

	// Output initial collagen fiber orientation
	fileout.writeCollFib(GetGrid(), true);

	return true;
}

//-----------------------------------------------------------------------------
// Initialize the nodal ECM values
bool FEAngio::InitECMDensity()
{
	int NN = m_grid.Nodes();
	vector<double> density(NN, 0.0);

	if (m_grid.m_coll_den == 0.0)
	{
		// get the material
		FEMaterial* pm = m_fem.GetMaterial(0);
		FEMaterial* pmat = pm->FindComponentByType("angio");
		if (pmat == 0) return false;

		if (CreateDensityMap(density, pmat) == false) return false;
	}
	else
	{
		for (int i=0; i<NN; ++i) density[i] = m_grid.m_coll_den;
	}

	// assign ECM density
	for (int i = 0; i < NN; ++i)								
	{
		Node& node = m_grid.GetNode(i);

		node.m_ecm_den0 = density[i];	
		node.m_ecm_den = node.m_ecm_den0;
	}

	return true;
}

//-----------------------------------------------------------------------------
// Initialize nodal collagen fiber directions
bool FEAngio::InitCollagenFibers()
{
	int NN = m_grid.Nodes();
	vector<vec3d> fiber;
	fiber.resize(NN, vec3d(0,0,0));

	switch (m_matrix_cond)
	{
	case 0: // random orientation
		{
			for (int i=0; i<NN; ++i)
			{
				fiber[i] = vrand();
			}
		}
		break;
	case 3:	// from element's local coordinates
		{
			// get the material
			FEMaterial* pm = m_fem.GetMaterial(0);
			FEMaterial* efd = pm->FindComponentByType("EFD neo-Hookean");
			if (efd == 0) return false;
			else
			{
				if (CreateFiberMap(fiber, efd) == false) return false;
			}
		}
	}

	// assign collagen fibers
	for (int i=0; i<NN; ++i)
	{
		Node& node = m_grid.GetNode(i);

		vec3d v = fiber[i];

		// flatten if requested
		if (m_bzfibflat == 1) v.z *= 0.25;

		// normalize the vector
		v.unit();

		// assign the node
		node.m_collfib0 = v;
		node.m_collfib = node.m_collfib0;
	}

	return true;
}

//-----------------------------------------------------------------------------
// This function runs the angiongenesis simulation code. 
// It basically repeats for each time step the following: first it does a growth step
// followed by an FE simulation (i.e. FEBio solve). 
bool FEAngio::Run()
{
	// Solve for initial step
	if (RunFEM() == false) return false;

	// Apply sprout forces at active growth tips
	apply_sprout_forces(1, 0.5);						// Apply sprout forces to the mesh
	adjust_mesh_stiffness();							// Adjust the stiffness of the mesh based on microvessel volume

	// This is the main time loop
	while (m_time.t < m_time.maxt)
	{	
		// Determine the current time step and update time
		updateTime();

		// Grow the culture (Elongation, Branching, and Anastomosis)
		m_pCult->Grow(m_time);
    	
		// Uncomment this if not using a composite material model.
		// This will update the stiffness of any element that contains microvessel segements 
//		adjust_mesh_stiffness();
		
		// Discretize the growth step over N steps in quasi-time to produce smoother results (default N = 2)
		if (Subgrowth(m_sub_cycles) == false) break;
		
		// Output time information	
		fileout.save_time(*this);
		
		// Print the status of angio3d to the user    
		fileout.printStatus(*this);
	}

	// generate all output
	Output();

	return true;
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{
	//// ANGIO3D - Generate output files
	output_params();										// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	
	//feangio.fileout.printsproutnodes(feangio.sprout_nodes);		// Output the sprout locations
	
	fileout.dataout(*this);								// Output data file
		
	fileout.writeCollFib(GetGrid(), false);				// Output final collagen fiber orientation

	fileout.writeECMDen(GetGrid());						// Output final matrix density
}

//-----------------------------------------------------------------------------
// This runs the FE model and updates the grid and other data that depends on
// the FE solution.
bool FEAngio::RunFEM()
{
	// increase the FE_state counter
	FE_state++;

	// Reset some parameters for FEBio
	FEAnalysis* pstep = m_fem.GetCurrentStep();
	pstep->m_dt         = pstep->m_dt0;
	pstep->m_ntimesteps = 0;
	pstep->m_iteopt     = 100;
	pstep->m_maxretries = 10;

	// Solve the FE problem
	if (m_fem.Solve() == false) return false;

	// Update the position of grid nodes using the solution from FEBio.
	// and dependent variables.
	UpdateGrid();

	// Save the current vessel state
	fileout.save_vessel_state(*this);

	return true;
}

//-----------------------------------------------------------------------------
// Step growth through a series of sub-steps to produce smoother model results.
bool FEAngio::Subgrowth(int sub_steps)
{
	update_sprout_stress_scaling();

	// Iterate through the number of substeps...
	for (int k = 1; k <= sub_steps; k++)
	{
		// Update the value of the subgrowth scaling factor
		double subgrowth_scale = ((double)k/(double)sub_steps);

		// do the sub-growht step
		m_pCult->SubGrowth(subgrowth_scale);

		// Update the positions of the body forces
		update_body_forces(1.0);

		// Run the FE analysis
		if (RunFEM() == false) return false;
	}
	
	// Remove bad segments
//	removeErrors();

	return true;
}

//-----------------------------------------------------------------------------
// Use the displacement field from the FE solution to update microvessels into the current configuration
void FEAngio::DisplaceVessels()
{
	// loop over all fragments
	for (SegIter it = m_pCult->SegmentBegin(); it != m_pCult->SegmentEnd(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		// Iterate through both segment tips
		for (int k=0; k<2; ++k)
		{
			// get the tip
			Segment::TIP& tip = it->tip(k);

			// Update position
			assert(tip.pt.nelem >= 0);
			tip.rt = m_grid.Position(tip.pt);
		}
		
		// Recalculate the segment's length and unit vector based on it's new position
		it->Update();
	}

	// Update the total vascular length within the simulation   
    UpdateTotalLength();
}   


/*
//-----------------------------------------------------------------------------
// Use the displacement field from the FE solution to update microvessels into the current configuration
void FEAngio::DisplaceVessels()
{
	double xix = 0.; double xiy = 0.; double xiz = 0.;			// Position in the element's natural coordinates		
	double shapeF[8] = {0.};									// Array containing the shape function values at the segment's position

	// loop over all fragments
	for (SegIter it = m_pCult->SegmentBegin(); it != m_pCult->SegmentEnd(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		// Iterate through both segment tips
		for (int k=0; k<2; ++k)
		{
			// get the tip
			Segment::TIP& tip = it->tip(k);

			// get the element number
			int elem_num = tip.pt.nelem;								// Find the element that contains the segment tip
			assert(elem_num >= 0);

			// get the element this tip is in
			Elem& elem = m_grid.ebin[elem_num];

			// Get position to the current segment tip
			vec3d rt_old = tip.rt;
		    
			// Convert the position to natural coordinates
			m_grid.natcoord(xix, xiy, xiz, rt_old.x, rt_old.y, rt_old.z, elem_num);
		    
			// Obtain the values of the shape functions at this position
			m_grid.shapefunctions(shapeF, xix, xiy, xiz);
		    
			// Calculate the displacement vector by interpolating nodal displacement to the segment tip
			vec3d disp;
			disp.x = shapeF[0]*(*elem.n1).u.x + shapeF[1]*(*elem.n2).u.x + shapeF[2]*(*elem.n3).u.x + shapeF[3]*(*elem.n4).u.x + shapeF[4]*(*elem.n5).u.x + shapeF[5]*(*elem.n6).u.x + shapeF[6]*(*elem.n7).u.x + shapeF[7]*(*elem.n8).u.x;
			disp.y = shapeF[0]*(*elem.n1).u.y + shapeF[1]*(*elem.n2).u.y + shapeF[2]*(*elem.n3).u.y + shapeF[3]*(*elem.n4).u.y + shapeF[4]*(*elem.n5).u.y + shapeF[5]*(*elem.n6).u.y + shapeF[6]*(*elem.n7).u.y + shapeF[7]*(*elem.n8).u.y;
			disp.z = shapeF[0]*(*elem.n1).u.z + shapeF[1]*(*elem.n2).u.z + shapeF[2]*(*elem.n3).u.z + shapeF[3]*(*elem.n4).u.z + shapeF[4]*(*elem.n5).u.z + shapeF[5]*(*elem.n6).u.z + shapeF[6]*(*elem.n7).u.z + shapeF[7]*(*elem.n8).u.z;
		    		    
			// Calculate the weighted displacement vector
			vec3d weighted_disp = disp*phi_stiff_factor;
		
			// Update the segment tip position using the new weighted displacement vector
			tip.rt += rt_old + weighted_disp;
			
			// If using the weighted displacement vector causes the segment to move outside the mesh...
			if (m_grid.findelem(tip.rt) == -1)
			{	
				// Update using the full displacement vector instead
				tip.rt = rt_old + disp;
			}
		}
		
		// Recalculate the segment's length and unit vector based on it's new position
		it->Update();
	}

	// Update the total vascular length within the simulation   
    UpdateTotalLength();
}   
*/

//-----------------------------------------------------------------------------
// Apply sprout forces to the mesh for each active vessel tip
void FEAngio::apply_sprout_forces(int load_curve, double scale)
{
	double magnitude = scale*m_sproutf;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

	//#pragma omp parallel for
	for (list<SegIter>::iterator tip_it = m_pCult->m_active_tips.begin(); tip_it != m_pCult->m_active_tips.end(); ++tip_it)		// For each active growth tip...
	{
		Segment& seg = (*(*tip_it));												// Obtain the growth tip

		if (seg.tip(0).bactive)
		{
			vec3d tip = seg.tip(0).rt;												// Obtain the position of the active tip
			
			// Calculate the directional unit vector of the sprout (notice negative sign)
			vec3d sprout_vect = -seg.uvect();

			(*tip_it)->tip(0).bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, load_curve);				// Create a new body force, set the tips body force ID
		}
		
		if (seg.tip(1).bactive)
		{	
			vec3d tip = seg.tip(1).rt;												// Obtain the position of the active tip
			
			// Calculate the directional unit vector of the sprout
			vec3d sprout_vect = seg.uvect();

			(*tip_it)->tip(1).bdyf_id = create_body_force(sprout_vect, tip.x, tip.y, tip.z, magnitude, 1.0/m_tip_range, load_curve);				// Create a new body force, set the tips body force ID
		}
	}

	return;
}

//-----------------------------------------------------------------------------
// Add a new body force entry into the body force field applyied to the mesh
int FEAngio::create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, int load_curve)
{
	total_bdyf++;							// Iterate the total body force counter							

	if (m_pmat)
	{
		m_pmat->AddSprout(vec3d(xpt, ypt, zpt), sprout_vect);
		return m_pmat->Sprouts() - 1;
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));					// Obtain the body force class
		pbf->AddSprout(vec3d(xpt, ypt, zpt),sprout_vect);												// Add a new component to the body force for this active sprout tip
		return pbf->Sprouts() - 1;																		// Return the ID number for the body force
	}
}


//-----------------------------------------------------------------------------
// Update the sprout forces after a deformation
void FEAngio::update_body_forces(double scale)
{
	vec3d sprout_vect;												// Sprout direction vector

	vec3d tip(0,0,0);
	double magnitude = scale*m_sproutf;								// Magnitude of the sprout force

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

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
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));			// Obtain the sprout body force field
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

	FEMesh& mesh = m_fem.GetMesh();									// Obtain the FE mesh

	//#pragma omp parallel for
	for (SegIter frag_it = m_pCult->SegmentBegin(); frag_it != m_pCult->SegmentEnd(); ++frag_it)		// Iterate through each segment in the model...
	{
		const Segment& seg = (*frag_it);								// Obtain the segment, keep it constant to prevent changes

		if (((seg.tip(0).bactive) || (seg.tip(0).BC == 1)) && (seg.tip(0).bdyf_id >= 0)){		  // Turn on the body force for any active -1 segment OR -1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[0] == -1) && (seg.bdyf_id[0] >= 0)){									// Turn on the body force for any active -1 segment
			tip = seg.tip(0).rt;																	// Obtain the tip position

			sprout_vect = seg.tip(0).rt - seg.tip(1).rt;												// Calculate the sprout directional vector
			sprout_vect = sprout_vect/sprout_vect.norm();			

			update_angio_sprout(seg.tip(0).bdyf_id, true, tip, sprout_vect);
			}
		
		if (((seg.tip(1).bactive) || (seg.tip(1).BC == 1)) && (seg.tip(1).bdyf_id >= 0)){		  // Turn on the body force for any active +1 segment OR +1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[1] == 1) && (seg.bdyf_id[1] >= 0)){									// Turn on the body force for any active +1 segment
			tip = seg.tip(1).rt;																	// Obtain the tip position
			
			sprout_vect = seg.tip(1).rt - seg.tip(0).rt;												// Calculate the sprout directional vector
			sprout_vect.unit();
			update_angio_sprout(seg.tip(1).bdyf_id, true, tip, sprout_vect);
			}
	}
	
	return;
}

//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Update the grid after a deformation using the FE mesh
void FEAngio::UpdateGrid()
{
	FEMesh& mesh = m_fem.GetMesh();

	// loop over all nodes
	int NN = m_grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		// Update the grid node to the current position of the FE mesh
		m_grid.GetNode(i).rt = mesh.Node(i).m_rt;
	}

	// Update the volume for each element in the grid
	m_grid.update_grid_volume();

	// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)
	update_ECM();

	// Displace the microvessels using the solution from FEBio
	DisplaceVessels();
}

//-----------------------------------------------------------------------------
// Update the ECM field after a deformation
void FEAngio::update_ECM()
{
	// Natural coordinates of each node within the element
	double LUT[8][3] = {
		{-1., -1. ,-1.},
		{ 1., -1. ,-1.},
		{-1.,  1. ,-1.},
		{ 1.,  1. ,-1.},
		{-1., -1. , 1.},
		{ 1., -1. , 1.},
		{-1.,  1. , 1.},
		{ 1.,  1. , 1.}
	};

	// reset nodal data
	int NN = m_grid.Nodes();
	for (int i=0; i<NN; ++i)
	{
		Node& ni = m_grid.GetNode(i);
		ni.m_ntag = 0;
		ni.m_collfib = vec3d(0,0,0);
		ni.m_ecm_den = 0.0;
	}

	// For each element within the grid...
	int NE = m_grid.Elems();
	for (int i = 0; i < NE; ++i)
	{
		// Obtain the element
		Elem& elem = m_grid.GetElement(i);
		
		// For each node in the element...
		for (int j=0; j<8; j++)
		{
			// get the node
			Node* nj = elem.GetNode(j);

			// get the ecm density and collagen fiber
			double ecm_den  = nj->m_ecm_den0;
			vec3d coll_fib = nj->m_collfib0;
			
			// Calculate the deformation gradient tensor and jacobian at the node
			mat3d F = m_grid.calculate_deform_tensor(elem, LUT[j][0], LUT[j][1], LUT[j][2]);
			double Jacob = F.det();
			assert(Jacob > 0.0);
			
			// Update the collagen fiber orientation vector into the current configuration using F		
			coll_fib = F*coll_fib;
			coll_fib.unit();

			// Update matrix density using the Jacobian
			ecm_den = ecm_den/Jacob;

			// accumulate fiber directions and densities
			nj->m_collfib += coll_fib;
			nj->m_ecm_den += ecm_den;

			// increment counter
			nj->m_ntag++;
		}
	}

	// normalize fiber vector and average ecm density
	for (int i = 0; i < NN; ++i)
	{
		Node& ni = m_grid.GetNode(i);
		assert(ni.m_ntag > 0);
		ni.m_ecm_den /= (double) ni.m_ntag;
		ni.m_collfib.unit();
	}
	
	//update_ecm_den_grad();										  // Update the ECM density gradient based on the solution from FEBio
		
	return;
}

///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////
// TODO: vessel lengths are always positive now, so we need to fix the logic here.
void FEAngio::adjust_mesh_stiffness()
{
	if (comp_mat == 0)													// If a composite consitutive model isn't being used, exit
		return;

	Grid& grid = m_grid;
		
	int elem_num = 0;													// Element number
	vec3d vess_vect;													// Vessel vector

	int NE = m_grid.Elems();
	for (int i = 0; i < NE; ++i)
	{							
		Elem& el = grid.GetElement(i);
		el.alpha = 0.;						// Set the vessel volume fraction, alpha, to zero
		el.fiber_orient = vec3d(0,0,0);		// Set the vessel orientation vector to 0 (this is the element fiber_orient vector, which does not contain collagen fiber orientation information but rather the material fiber direction for a transversely isotropic material model)
	}
	
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1/(double)Nsub;									// Calculate the subdivision scale 
	vec3d mid;
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

	for (SegIter frag_it = m_pCult->SegmentBegin(); frag_it != m_pCult->SegmentEnd(); ++frag_it)		// For each segment...
	{
		Segment subunit;												// Segment subdivision placeholder
	
		Segment& seg = (*frag_it);												// Obtain the segment
		
		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length() > 0.){											// If it's a +1 segment...
					subunit.tip(0).rt = seg.tip(0).rt;
				}
				else{															// If it's a -1 segment...
					subunit.tip(1).rt = seg.tip(1).rt;
//					subunit.m_length = -1.;
				}
		}
			
			// Calculate the subdivision
			if (seg.length() > 0.){										// If it's a +1 segment...			
				subunit.tip(1).rt = subunit.tip(0).rt + seg.uvect()*(sub_scale*seg.length());     
			}
			else{														// If it's a -1 segment...
				subunit.tip(0).rt = subunit.tip(1).rt + seg.uvect()*(sub_scale*seg.length());     
			}

			subunit.Update();										// Find the length of the subdivision
			
			mid = (subunit.tip(1).rt + subunit.tip(0).rt)*0.5;

			elem_num = m_grid.findelem(mid);				// Find the element that the midpoint is within

			// Calculate the orientation of the subdivision
			if (seg.length() > 0.){										// If it's a +1 segment...
				vess_vect = subunit.tip(1).rt - subunit.tip(0).rt;
			}
			else{														// If it's a -1 segment...
				vess_vect = subunit.tip(0).rt - subunit.tip(1).rt;
			}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
				vess_vect = vess_vect/vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				Elem& el = m_grid.GetElement(elem_num);

				elem_volume = el.volume;					// Calculate the volume of the element

				subunit_volume = pi*(m_vessel_width/2.)*(m_vessel_width/2.)*fabs(subunit.length());		// Find the volume of the subdivision
				volume_fraction = subunit_volume/elem_volume;				// Calculate the volume fraction

				el.alpha = el.alpha + volume_fraction;	// Add the volume fraction for each subdivision to alpha
			
				// Calculate the vessel orientation vector 
				if ((el.fiber_orient.x == 0) && (el.fiber_orient.y == 0) && (el.fiber_orient.z == 0)){	// If the vessel orientation vector hasn't been assigned yet...
					el.fiber_orient = vess_vect;			// Set the vessel orientation vector					
					el.fiber_orient.z = vess_vect.z;
				}
				else{														// If it has been...	
					el.fiber_orient.x = (el.fiber_orient.x + vess_vect.x)/2;	// Average together the vessel orientation vector
					el.fiber_orient.y = (el.fiber_orient.y + vess_vect.y)/2;
					el.fiber_orient.z = (el.fiber_orient.z + vess_vect.z)/2;}
			}
			
			// Set the origin of the next subdivision to the end of the current one
			if (seg.length() > 0.){
				subunit.tip(0).rt = subunit.tip(1).rt;
			}
			else{
				subunit.tip(1).rt = subunit.tip(0).rt;
			}
		}
	}
	
	double alpha = 0.;													// Volume fraction for the composite material model
	vec3d e1; vec3d e2; vec3d e3;										// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
	
	FEMesh& mesh = m_fem.GetMesh();										// Get the FE mesh
	int J = mesh.Domains();												// Find the number of domains within the mesh
    int num_elem = 0;

	for (int k = 0; k < J; ++k){
		FEDomain& d = mesh.Domain(k);										// Obtain the domain
		
		for (int j = 0; j < d.Elements(); ++j)								// For each element within the domain...
		{
			FEElement& e = d.ElementRef(j);										// Obtain the element from the domain
			int nint = e.GaussPoints();											// Obtain the number of gauss points
		
			Elem& eg = grid.GetElement(num_elem);
			alpha = eg.alpha;											// Obtain alpha from the grid element

			// Set e1 to the vessel orientation vector
			e1 = eg.fiber_orient;

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

//-----------------------------------------------------------------------------
// Calculate the density gradient for each element (this function may not work properly)
void FEAngio::update_ecm_den_grad()
{
	assert(false);
/*	vec3d elem_den_grad;
	double ex = 0.; double ey = 0.; double ez = 0.;	
	vec3d dN1; vec3d dN2; vec3d dN3; vec3d dN4; vec3d dN5; vec3d dN6; vec3d dN7; vec3d dN8;
	
	dN1 = m_grid.shapefun_d1(ex, ey, ez, 1);
    dN2 = m_grid.shapefun_d1(ex, ey, ez, 2);
	dN3 = m_grid.shapefun_d1(ex, ey, ez, 3);
	dN4 = m_grid.shapefun_d1(ex, ey, ez, 4);
	dN5 = m_grid.shapefun_d1(ex, ey, ez, 5);
	dN6 = m_grid.shapefun_d1(ex, ey, ez, 6);
	dN7 = m_grid.shapefun_d1(ex, ey, ez, 7);
	dN8 = m_grid.shapefun_d1(ex, ey, ez, 8);
	
	Grid& grid = m_grid;
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
*/
}


void FEAngio::update_sprout_stress_scaling()
{
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	if (m_pmat)
		m_pmat->scale = y0 + a/(1 + exp(-(m_time.t - x0)/b));
	
	return;
}

///////////////////////////////////////////////////////////////////////
// enfore_fiber_BCS
///////////////////////////////////////////////////////////////////////

void FEAngio::enforce_fiber_BCS(Node &node, bool circ)
{
	Grid& grid = m_grid;
	if (circ == true){
		vec3d r;
		r.x = node.rt.x - grid.xrange[1];
		r.y = node.rt.y - grid.yrange[1];

		r = r/r.norm();

		mat3d I(1.0);
		mat3d M = r&r;

		vec3d r_new; r_new = (I - M)*node.m_collfib; r_new = r_new/r_new.norm();

		node.m_collfib.x = r_new.x; node.m_collfib.y = r_new.y; node.m_collfib.z = r_new.z;

		return;
	}
		
	if (m_cgelbc == 'l'){
		//Node on the front face
		if (node.rt.y == grid.yrange[0])
			node.m_collfib.y = 0;

		// Node on the back face
		//if (node.rt.y == grid.yrange[1])
			//node.collfib.y = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.m_collfib.z = 0;

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
			node.m_collfib.x = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.m_collfib.z = 0;

		// Node on the top face
		//if (node.rt.z == grid.zrange[1])
			//node.collfib.z = 0;
	}

	if (m_cgelbc == 'u'){
		//Node on the front face
		if (node.rt.y == grid.yrange[0])
			node.m_collfib.y = 0;

		// Node on the right face
		//if (node.rt.x == grid.xrange[0])
			//node.collfib.x = 0;

		// Node on the back face
		//if (node.rt.y == grid.yrange[1])
			//node.collfib.y = 0;

		// Node on the left face
		if (node.rt.x == grid.xrange[0])
			node.m_collfib.x = 0;

		// Node on the bottom face
		if (node.rt.z == grid.zrange[0])
			node.m_collfib.z = 0;

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

	// get the mesh
	FEMesh& mesh = pmat->GetFEModel()->GetMesh();

	// initialize the fiber array
	int N = mesh.Nodes();
	fiber.assign(N, vec3d(0,0,0));

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

//-----------------------------------------------------------------------------
// Update the time step size and current time value.
void FEAngio::updateTime()
{
	double dt = m_dtA*exp(m_dtB/(m_ntime + m_dtC));
    m_ntime += 1;

    if (dt > (m_time.maxt - m_time.t) && (m_time.maxt - m_time.t) > 0)       // If dt is bigger than the amount of time left...
		dt = m_time.maxt - m_time.t;                                       // then just set dt equal to the amount of time left (maxt-t)

	// Update time
	m_time.dt = dt;
    m_time.t += dt; 
}

//-----------------------------------------------------------------------------
// Calculates and stores the total length of all vessels.
void FEAngio::UpdateTotalLength()
{
    m_total_length = 0.;
    for (SegIter it = m_pCult->SegmentBegin(); it != m_pCult->SegmentEnd(); ++it)
    {
        m_total_length += it->length();
    }
}

///////////////////////////////////////////////////////////////////////
// removeErrors
///////////////////////////////////////////////////////////////////////

void FEAngio::removeErrors()
{
	if (kill_off == false){
		double length_limit = m_pCult->m_d;
    
		for (SegIter it = m_pCult->SegmentBegin(); it != m_pCult->SegmentEnd(); ++it){
			if (fabs(it->length()) >= 2.0*length_limit){
				it->death_label = -7;
				vec3d& r0 = it->tip(0).rt;
				vec3d& r1 = it->tip(1).rt;
//				fprintf(killed_segs_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->death_label);
//				it = m_pCult->m_frag.erase(it);
				assert(false);
			}
			if (it == m_pCult->SegmentEnd())
				return;
		}
                
		return;
	}
}
