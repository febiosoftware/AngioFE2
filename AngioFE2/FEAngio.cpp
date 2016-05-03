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
#include <FECore/FEAnalysis.h>
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

	// boundary conditions
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
FEAngioMaterial* FEAngio::FindAngioMaterial(FEMaterial* pm)
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

	// create the grid based on the FEBio mesh
	if (m_grid.Init() == false) return false;

	// assign ECM densities to grid nodes
	if (InitECMDensity() == false) return false;

	// assign collagen fibers to grid nodes
	if (InitCollagenFibers() == false) return false;

	// initialize culture
	// NOTE: must be done after InitECMDensity() and InitCollagenFibers().
	// Initialize Culture class
	if (m_pCult->Init() == false) return false;

	// Init all the FE stuff
	if (InitFEM() == false) return false;
	
	// start timer
	time(&m_start);

	return true;
}

//-----------------------------------------------------------------------------
double FEAngio::RunTime()
{
	time_t stop;
    time(&stop);
	return (double) difftime(stop, m_start);
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

	// register the callback
	m_fem.AddCallback(FEAngio::feangio_callback, CB_UPDATE_TIME | CB_MAJOR_ITERS | CB_SOLVED, this);

	// Do the model initialization
	if (m_fem.Init() == false) return false;

	// apply the intial sprouts
	CreateSprouts(0.5);

	// Adjust the stiffness of the mesh based on microvessel volume
	adjust_mesh_stiffness();

	// only output to the logfile (not to the screen)
	felog.SetMode(Logfile::FILE_ONLY);

	// --- Output initial state of model ---

	// Output initial microvessel state
	fileout.save_vessel_state(*this);

	// save active tips
	fileout.save_active_tips(*this);

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
void FEAngio::OnCallback(FEModel* pfem, unsigned int nwhen)
{
	FEModel& fem = *pfem;

	if (nwhen == CB_UPDATE_TIME)
	{
		// grab the time information
		m_time.t = fem.m_ftime;
		m_time.dt = fem.GetCurrentStep()->m_dt;

		// do a growth step
		m_pCult->Grow(m_time);

		// update sprout stress scaling
		update_sprout_stress_scaling();

		// Update the positions of the body forces
		UpdateSprouts(1.0);
	}
	else if (nwhen == CB_MAJOR_ITERS)
	{
		// update the grid data
		m_grid.Update();

		// update the culture
		m_pCult->Update();

		++FE_state;

		// Save the current vessel state
		fileout.save_vessel_state(*this);

		// save active tips
		fileout.save_active_tips(*this);

		// Output time information	
		fileout.save_time(*this);
		
		// Print the status of angio3d to the user    
		fileout.printStatus(*this);
	}
	else if (nwhen == CB_SOLVED)
	{
		// do the final output
		Output();
	}
}

//-----------------------------------------------------------------------------
// Generate output. This is called at the end of Run().
// Note that some output is not written here. E.g. the vessel state is written
// after each successful FE run.
void FEAngio::Output()
{
	// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	fileout.output_params(*this);
	
	// Output data file
	fileout.dataout(*this);
		
	// Output final collagen fiber orientation
	fileout.writeCollFib(GetGrid(), false);

	// Output final matrix density
	fileout.writeECMDen(GetGrid());
}

//-----------------------------------------------------------------------------
int FEAngio::Sprouts()
{
	if (m_pmat) return m_pmat->Sprouts();
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));
		return pbf->Sprouts();
	}
}

//-----------------------------------------------------------------------------
// Apply sprout forces to the mesh for each active vessel tip
void FEAngio::CreateSprouts(double scale)
{
	double magnitude = scale*m_sproutf;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_pCult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)
	{
		Segment::TIP& tip = *(*tip_it);
		if (tip.bactive)
		{
			// get the tip
			const vec3d& r = tip.pos();
			
			// get the directional unit vector of the tip
			const vec3d& u = tip.u;

			if (m_pmat)
			{
				m_pmat->AddSprout(r, u);
				tip.bdyf_id = m_pmat->Sprouts() - 1;
			}
			else
			{
				FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));
				pbf->AddSprout(r, u);
				tip.bdyf_id = pbf->Sprouts() - 1;
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Update the sprout forces after a deformation
void FEAngio::UpdateSprouts(double scale)
{
	double magnitude = scale*m_sproutf;								// Magnitude of the sprout force
	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (m_time.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (m_time.t < 4.0)
		magnitude = (1.0/4.0)*m_time.t*scale;

	if (m_pmat)
	{
		m_pmat->ClearSprouts();
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(m_fem.GetBodyLoad(0));			// Obtain the sprout body force field
		pbf->ClearSprouts();
		FEParameterList& pl = pbf->GetParameterList();										// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);												// Get the sprout force magnitude parameter
		pa->value<double>() = magnitude*m_sproutf;													// Set the sprout force magnitude parameter
	}

	// Obtain the FE mesh
	FEMesh& mesh = m_fem.GetMesh();

	//#pragma omp parallel for
	const SegmentTipList& tip_list = m_pCult->GetActiveTipList();
	for (ConstTipIter tip_it = tip_list.begin(); tip_it != tip_list.end(); ++tip_it)		// Iterate through each segment in the model...
	{
		const Segment::TIP& tip= *(*tip_it);
		assert(tip.bactive);
		assert(tip.bdyf_id >= 0);

		// TODO: What to do with BC==1? Currently, tips that stop growing after hitting boundary
		//       are no longer active. We should still add a sprout for those
		if      (m_pmat) m_pmat->AddSprout(tip.pos(), tip.u);
		else if (m_pbf ) m_pbf->AddSprout (tip.pos(), tip.u);
	}
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

	const SegmentList& seg_list = m_pCult->GetSegmentList();
	for (ConstSegIter frag_it = seg_list.begin(); frag_it != seg_list.end(); ++frag_it)		// For each segment...
	{
		Segment subunit;												// Segment subdivision placeholder
	
		const Segment& seg = (*frag_it);												// Obtain the segment
		
		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length() > 0.){											// If it's a +1 segment...
					subunit.tip(0).pt.r = seg.tip(0).pos();
				}
				else{															// If it's a -1 segment...
					subunit.tip(1).pt.r = seg.tip(1).pos();
//					subunit.m_length = -1.;
				}
		}
			
			// Calculate the subdivision
			vec3d v = seg.uvect();
			if (seg.length() > 0.){
				subunit.tip(1).pt.r = subunit.tip(0).pos() + v*(sub_scale*seg.length());     
			}
			else{														// If it's a -1 segment...
				subunit.tip(0).pt.r = subunit.tip(1).pos() + v*(sub_scale*seg.length());     
			}

			subunit.Update();										// Find the length of the subdivision
			
			mid = (subunit.tip(1).pos() + subunit.tip(0).pos())*0.5;

			elem_num = m_grid.findelem(mid);				// Find the element that the midpoint is within

			// Calculate the orientation of the subdivision
			if (seg.length() > 0.){										// If it's a +1 segment...
				vess_vect = subunit.tip(1).pos() - subunit.tip(0).pos();
			}
			else{														// If it's a -1 segment...
				vess_vect = subunit.tip(0).pos() - subunit.tip(1).pos();
			}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
				vess_vect = vess_vect/vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				Elem& el = m_grid.GetElement(elem_num);

				elem_volume = el.m_volume;					// Calculate the volume of the element

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
				subunit.tip(0).pt.r = subunit.tip(1).pos();
			}
			else{
				subunit.tip(1).pt.r = subunit.tip(0).pos();
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
				FEMaterialPoint& mp = *(e.GetMaterialPoint(n));
				FEAngioMaterialPoint* pt = FEAngioMaterialPoint::FindAngioMaterialPoint(&mp); // get the mixture material point
				pt->vessel_weight = alpha;
				pt->matrix_weight = 1.0 - alpha;

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

//-----------------------------------------------------------------------------
void FEAngio::update_sprout_stress_scaling()
{
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	if (m_pmat)
		m_pmat->scale = y0 + a/(1 + exp(-(m_time.t - x0)/b));
	
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
				FEAngioMaterialPoint* apoint = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
				FEElasticMaterialPoint& pt = *apoint->matPt->ExtractData<FEElasticMaterialPoint>();
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
				FEAngioMaterialPoint* angioPt = FEAngioMaterialPoint::FindAngioMaterialPoint(mpoint);
				den[n] = angioPt->m_D;
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
