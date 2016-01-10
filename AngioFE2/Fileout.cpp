///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "Culture.h"
#include "Data.h"
#include "angio3d.h"
#include "Grid.h"
#include "Elem.h"
#include "FEAngio.h"
#include <iostream>

#include "FEAngio.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Fileout::Fileout()
{
    logstream.open("out_log.ang");
    //stream3 = fopen("tracking.ang","wt");   // tracking.ang: time step, model time, total length in culture, number of branches in culture
	stream = fopen("out_data.ang","wt");                                        // data.ang: Store 3D coordinates of begining and end of each vessel segment
																				// as well as total length of the segment
	num_fe_timesteps = 48;
}

Fileout::~Fileout()
{
    logstream.close();
}

///////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// timestart
///////////////////////////////////////////////////////////////////////

void Fileout::timestart()
{
    time(&start);                                               // Start the timer
    return;
}



///////////////////////////////////////////////////////////////////////
// writeTracking
///////////////////////////////////////////////////////////////////////

void Fileout::writeTracking(Data &data)
{
    fprintf(stream3,"%-12.7f %-12.7f %-12.7f %-5i\n",data.dt,data.t,data.total_length,data.num_branches);   // Write to tracking.ang
    
    return;
}



///////////////////////////////////////////////////////////////////////
// closeTracking
///////////////////////////////////////////////////////////////////////

void Fileout::closeTracking()
{
    fclose(stream3);                                                        // Close stream to 'tracking.ang' (stream3) 
    
    return;
}



///////////////////////////////////////////////////////////////////////
// printStatus
///////////////////////////////////////////////////////////////////////

void Fileout::printStatus(Data &data)
{
    cout << endl << "Time: " << data.t << endl;                             // Print out current time to user
	//cout << "dt: " << data.dt << endl;
    cout << "Segments: " << data.nsegs << endl;                             // Print out current number of segments to user
	cout << "Total Length: " << data.total_length << endl;                  // Print out the current total length to user
	cout << "Branch Points: " << data.num_branches << endl;                 // Print out the current number of branches to user
	cout << "Anastomoses: " << data.num_anastom << endl << endl;            // Print out the current number of anastomoses to user
    
    logstream << endl << "Time: " << data.t << endl;                        // Print out current time to log file
	//logstream << "dt: " << data.dt << endl;
    logstream << "Segments: " << data.nsegs << endl;                        // Print out current number of segments to log file
	logstream << "Total Length: " << data.total_length << endl;             // Print out the current total length to log file
	logstream << "Branch Points: " << data.num_branches << endl;            // Print out the current number of branches to log file
	logstream << "Anastomoses: " << data.num_anastom << endl << endl;       // Print out the current number of anastomoses to log file
        
    return;
}



///////////////////////////////////////////////////////////////////////
// dataout
///////////////////////////////////////////////////////////////////////

void Fileout::dataout(FEAngio &feangio)
{
    writeData(feangio);                                    // Create and write to 'data.ang'
	
    //writeGrid(data, grid);                              // Create and write to 'grid.ang'
   
    //writeNodes(angio.data, angio.grid);
    
    //writeEconn(angio.data, angio.grid);
    
    //writeBC(angio.grid);
    
    //writeGrad(data, grid);                              // Create and write to 'grad.ang'
    
    //writeAngle(frag);                                   // Create and write to 'angle.ang'

    printtime();                                        // Display the run-time to the user

    return;
}



///////////////////////////////////////////////////////////////////////
// writeData
///////////////////////////////////////////////////////////////////////

void Fileout::writeData(FEAngio &feangio)
{
	list<Segment>::iterator it;
	
	fprintf(stream,"%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-5s %-5s\n","Time","X1","Y1","Z1","X2","Y2","Z2","Length","Vess","Label");  // Write column labels to data.ang
	
	for (it = feangio.cult.m_frag.begin(); it != feangio.cult.m_frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		vec3d& r0 = it->m_tip[0].rt;
		vec3d& r1 = it->m_tip[1].rt;
		fprintf(stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,r0.x,r0.y,r0.z,r1.x,r1.y,r1.z,it->length,it->seg_num,it->label);
	}
	fclose(stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}


///////////////////////////////////////////////////////////////////////
// writeNodes
///////////////////////////////////////////////////////////////////////

void Fileout::writeNodes(Data &data, Grid &grid)
{
    /// File output: 'nodes.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_nodes.ang","wt");                                       
	Node node;

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i){
	    node = grid.nodes[i];
	    fprintf(stream2, "%-5.2i %-12.7f %-12.7f %-12.7f\n", node.id, node.rt.x, node.rt.y, node.rt.z);
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;
}



///////////////////////////////////////////////////////////////////////
// writeEconn
///////////////////////////////////////////////////////////////////////

void Fileout::writeEconn(Data &data, Grid &grid)
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
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].collfib.x, grid.nodes[i].collfib.y, grid.nodes[i].collfib.z);
	}

	fclose(node_stream);

	return;
}



void Fileout::writeECMDen(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den.ang","wt");

	//fprintf(node_stream,"%-5s %-12s %-12s %-12s %-12s %-12s\n","ID"," X "," Y "," Z ","ECM_DEN","ECM_DEN0");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].ecm_den, grid.nodes[i].ecm_den0);
	}

	fclose(node_stream);

	return;
}



void Fileout::writeECMDenGrad(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den_grad.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.nodes[i].id, grid.nodes[i].rt.x, grid.nodes[i].rt.y, grid.nodes[i].rt.z, grid.nodes[i].ecm_den_grad.x, grid.nodes[i].ecm_den_grad.y, grid.nodes[i].ecm_den_grad.z);
	}

	fclose(node_stream);

	return;
}

void Fileout::writeECMDenStore(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_density_store.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i ", grid.nodes[i].id+1);

		for (int j = 0; j < grid.nodes[i].ecm_den_store.size(); j++)
			fprintf(node_stream,"%-12.7f ",grid.nodes[i].ecm_den_store[j]);
		
		fprintf(node_stream,"\n");
	}

	fclose(node_stream);

	return;
}

void Fileout::writeECMFibrilStore(Grid &grid)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_fibril_store.ang","wt");	

	int NN = grid.Nodes();
	for (int i = 0; i < NN; ++i)
	{
		fprintf(node_stream,"%-5.2i ", grid.nodes[i].id+1);

		for (int j = 0; j < grid.nodes[i].ecm_den_store.size(); j++)
			fprintf(node_stream,"%-12.7f %-12.7f %-12.7f ",grid.nodes[i].ecm_fibril_store[j].x,grid.nodes[i].ecm_fibril_store[j].y,grid.nodes[i].ecm_fibril_store[j].z);
		
		fprintf(node_stream,"\n");
	}

	fclose(node_stream);

	return;
}

///////////////////////////////////////////////////////////////////////
// writeBC
///////////////////////////////////////////////////////////////////////

void Fileout::writeBC(Grid &grid)
{
    /// File output: 'eBC.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_eBC.ang","wt");                                       
	
	int BC_violate[6] = {0};

	int NE = grid.Elems();
	for (int i = 0; i < NE; ++i){
	    Elem& elem = grid.ebin[i];
	    
	    for (int j = 0; j < 6; ++j)
	        BC_violate[j] = 0;
	    
	    if ((elem.f1.BC == true) || (elem.f2.BC == true) || (elem.f3.BC == true) || (elem.f4.BC == true) || (elem.f5.BC == true) || (elem.f6.BC == true)){
	        if (elem.f1.BC == true)
	            BC_violate[0] = 1;
	        if (elem.f2.BC == true)
	            BC_violate[1] = 1;
	        if (elem.f3.BC == true)
	            BC_violate[2] = 1;
	        if (elem.f4.BC == true)
	            BC_violate[3] = 1;
	        if (elem.f5.BC == true)
	            BC_violate[4] = 1;
	        if (elem.f6.BC == true)
	            BC_violate[5] = 1;   
	        
	        fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, BC_violate[0], BC_violate[1], BC_violate[2], BC_violate[3], BC_violate[4], BC_violate[5]);}
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;

}


///////////////////////////////////////////////////////////////////////
// printtime
///////////////////////////////////////////////////////////////////////

void Fileout::printtime()
{
    time(&stop);                                                // Stop the timer
	
	t_seconds = (double) difftime(stop, start);                 // Calculate the simulation time in seconds
	
	cout << endl << "Simulation time: " << t_seconds << " seconds (" 
	    << floor(t_seconds/60) << " minutes)." << endl << endl;                // Show the user how long the simulation took (in seconds)
    
    logstream << endl << "Simulation time: " << t_seconds << " seconds (" << floor(t_seconds/60) << " minutes)." << endl << endl;  
	    
    return;
}

void Fileout::printrandseed(int randseed)
{
	logstream << endl << "Rand seed:" << randseed << endl << endl;
	return;
}

void Fileout::writeSegConn(list<Segment> &frag)
{
	list<Segment>::iterator it;
	
	///// File output: 'out_seg_conn.ang' /////
	
	FILE *stream;                                                           // Open stream to 'out_seg_conn.ang' (stream)
	stream = fopen("out_seg_conn.ang","wt");                                        
	
	fprintf(stream,"%-5s %-5s %-5s %-5s %-5s\n","SNum","Node1","     ","Node2","     ");  // Write column labels to data.ang

	for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		fprintf(stream,"%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i\n",it->seg_num,it->seg_conn[0][0],it->seg_conn[0][1],it->seg_conn[1][0],it->seg_conn[1][1]);  // Write to seg_conn.ang
	
	}
	
	fclose(stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}