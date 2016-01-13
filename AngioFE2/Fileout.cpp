///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "Culture.h"
#include "angio3d.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/FEModel.h>

#include "FEAngio.h"

using namespace std;

//-----------------------------------------------------------------------------
Fileout::Fileout()
{
    logstream.open("out_log.ang");
    //stream3 = fopen("tracking.ang","wt");   // tracking.ang: time step, model time, total length in culture, number of branches in culture
	m_stream = fopen("out_data.ang","wt");                                        // data.ang: Store 3D coordinates of begining and end of each vessel segment
																			// as well as total length of the segment}
	m_stream2 = fopen("out_vess_state.ang","wt");						// Open the stream for the vessel state data file		
	bf_stream = fopen("out_bf_state.ang","wt");						// Open the stream for the body force state data file

	time_stream = fopen("out_time.ang","wt");						// Open the stream for the time and state data file
	time_write_headers = true;										// Set the time and state data file to write the headers on its first output
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	fclose(m_stream2);
	fclose(bf_stream);
	fclose(time_stream);
}

//-----------------------------------------------------------------------------
void Fileout::writeTracking(FEAngio& angio)
{
	Culture& cult = angio.GetCulture();
    fprintf(m_stream3,"%-12.7f %-12.7f %-12.7f %-5i\n",angio.m_time.dt,angio.m_time.t,angio.m_total_length,cult.m_num_branches);   // Write to tracking.ang
    
    return;
}

//-----------------------------------------------------------------------------
void Fileout::closeTracking()
{
    fclose(m_stream3);                                                        // Close stream to 'tracking.ang' (stream3) 
    
    return;
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio)
{
	Culture& cult = angio.GetCulture();
    cout << endl << "Time: " << angio.m_time.t << endl;                             // Print out current time to user
	//cout << "dt: " << data.dt << endl;
    cout << "Segments: " << cult.Segments() << endl;                             // Print out current number of segments to user
	cout << "Total Length: " << angio.m_total_length << endl;                  // Print out the current total length to user
	cout << "Branch Points: " << cult.m_num_branches << endl;                 // Print out the current number of branches to user
	cout << "Anastomoses: " << cult.m_num_anastom << endl << endl;            // Print out the current number of anastomoses to user
    
    logstream << endl << "Time: " << angio.m_time.t << endl;                        // Print out current time to log file
	//logstream << "dt: " << data.dt << endl;
    logstream << "Segments: " << cult.Segments() << endl;                        // Print out current number of segments to log file
	logstream << "Total Length: " << angio.m_total_length << endl;             // Print out the current total length to log file
	logstream << "Branch Points: " << cult.m_num_branches << endl;            // Print out the current number of branches to log file
	logstream << "Anastomoses: " << cult.m_num_anastom << endl << endl;       // Print out the current number of anastomoses to log file
        
    return;
}

//-----------------------------------------------------------------------------
void Fileout::dataout(FEAngio &feangio)
{
    writeData(feangio);                                    // Create and write to 'data.ang'
	
    //writeGrid(data, grid);                              // Create and write to 'grid.ang'
   
    //writeNodes(angio.data, angio.grid);
    
    //writeEconn(angio.data, angio.grid);
    
    //writeBC(angio.grid);
    
    //writeGrad(data, grid);                              // Create and write to 'grad.ang'
    
    //writeAngle(frag);                                   // Create and write to 'angle.ang'

    printtime(feangio);                                        // Display the run-time to the user

    return;
}

//-----------------------------------------------------------------------------
void Fileout::writeData(FEAngio &angio)
{
	fprintf(m_stream,"%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-5s %-5s\n","Time","X1","Y1","Z1","X2","Y2","Z2","Length","Vess","Label");  // Write column labels to data.ang
	
	Culture& cult = angio.GetCulture();
	SegIter it;
	for (it = cult.SegmentBegin(); it != cult.SegmentEnd(); ++it)
	{
		vec3d& r0 = it->tip(0).rt;
		vec3d& r1 = it->tip(1).rt;
		fprintf(m_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->GetTimeOfBirth(),r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length(),it->m_nid,it->m_nseed);
	}
	fclose(m_stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeNodes(FEAngio& angio)
{
    Grid& grid = angio.GetGrid();
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_nodes.ang","wt");                                       
	Node node;

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i){
	    Node& node = grid.GetNode(i);
	    fprintf(stream2, "%-5.2i %-12.7f %-12.7f %-12.7f\n", node.m_id, node.rt.x, node.rt.y, node.rt.z);
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeEconn(FEAngio& angio)
{
/*    /// File output: 'econn.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_econn.ang","wt");                                       
	Elem elem;
	
	for (int i = 0; i < grid.Ne; ++i){
	    elem = grid.ebin[i];
	    fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, elem.n1.id, elem.n2.id, elem.n3.id, elem.n4.id, elem.n5.id, elem.n6.id,  elem.n7.id,  elem.n8.id);
	}  
	                                                                      	
	fclose(stream2); */                                                       
	
	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeCollFib(Grid &grid, bool initial)
{
	FILE *node_stream;
	
	if (initial == true)
		node_stream = fopen("out_coll_fib_init.ang","wt");
	else
		node_stream = fopen("out_coll_fib.ang","wt");

	//fprintf(node_stream,"%-5s %-12s %-12s %-12s %-12s %-12s\n","ID"," X "," Y "," Z ","THETA","ETA");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		Node& ni = grid.GetNode(i);
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", ni.m_id, ni.rt.x, ni.rt.y, ni.rt.z, ni.m_collfib.x, ni.m_collfib.y, ni.m_collfib.z);
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::writeECMDen(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den.ang","wt");

	//fprintf(node_stream,"%-5s %-12s %-12s %-12s %-12s %-12s\n","ID"," X "," Y "," Z ","ECM_DEN","ECM_DEN0");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		Node& ni = grid.GetNode(i);
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", ni.m_id, ni.rt.x, ni.rt.y, ni.rt.z, ni.m_ecm_den, ni.m_ecm_den0);
	}

	fclose(node_stream);

	return;
}

//-----------------------------------------------------------------------------
void Fileout::printtime(FEAngio& angio)
{
	time_t stop;
    time(&stop);                                                // Stop the timer
	
	double t_seconds = (double) difftime(stop, angio.m_start);                 // Calculate the simulation time in seconds
	
	cout << endl << "Simulation time: " << t_seconds << " seconds (" 
	    << floor(t_seconds/60) << " minutes)." << endl << endl;                // Show the user how long the simulation took (in seconds)
    
    logstream << endl << "Simulation time: " << t_seconds << " seconds (" << floor(t_seconds/60) << " minutes)." << endl << endl;  
}

//-----------------------------------------------------------------------------
void Fileout::printrandseed(int randseed)
{
	logstream << endl << "Rand seed:" << randseed << endl << endl;
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	fprintf(m_stream2,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length");  // Write column labels to out_vess_state.ang
	
	Culture& cult = angio.GetCulture();
	for (SegIter it = cult.SegmentBegin(); it != cult.SegmentEnd(); ++it)	// Iterate through all segments in frag list container (it)
	{
		vec3d& r0 = it->tip(0).rt;
		vec3d& r1 = it->tip(1).rt;
		fprintf(m_stream2,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",angio.FE_state,it->GetTimeOfBirth(),r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length());  // Write to out_vess_state.ang
	}
	
	return;
}

//-----------------------------------------------------------------------------
// Save positions of the body forces at the current time step (This function needs to be re-written)
void Fileout::save_bdy_forces(FEAngio& angio)
{
	FEModel& fem = angio.GetFEModel();
	int NBF = fem.BodyLoads();

	fprintf(bf_stream,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length"); 

	for (int i = 0; i < NBF; ++i)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
		FEParameterList& pl = pbf->GetParameterList();
		FEParam* pa = pl.Find("a");
		FEParam* prc = pl.Find("rc");

		if (pa && prc)
		{
			if (pa->value<double>() != 0.0)
				fprintf(bf_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",angio.FE_state, angio.m_time.t, prc->value<vec3d>().x, prc->value<vec3d>().y, prc->value<vec3d>().z, prc->value<vec3d>().x + 1.0, prc->value<vec3d>().y + 1.0, prc->value<vec3d>().z + 1.0, 1.73); 
		}
	}

	return;
}

//-----------------------------------------------------------------------------
// Save the current time information			
void Fileout::save_time(FEAngio& angio)
{
	if (time_write_headers == true){												// If this is the first time writing to out_time.ang
		fprintf(time_stream,"%-5s %-12s %-12s\n","State","Vess Time","FE Time");		// Print the column labels
		time_write_headers = false;}													// Turn off the headers flag
	
//	fprintf(time_stream,"%-5.2i %-12.7f %-12.7f\n",angio.FE_state, angio.m_time.t, ((double)angio.FE_state - 1.0)*angio.FE_time_step);	// Print out the FE state, the vessel growth model time, and the FE time
	fprintf(time_stream,"%-5.2i %-12.7f \n",angio.FE_state, angio.m_time.t);	// Print out the FE state, the vessel growth model time, and the FE time
	
	return;
}
