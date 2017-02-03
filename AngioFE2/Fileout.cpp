///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "Culture.h"
#include "Elem.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/FEModel.h>
#include "FEAngioMaterial.h"


using namespace std;

//-----------------------------------------------------------------------------
Fileout::Fileout()
{
    logstream.open("out_log.csv");
	//write the headers
	logstream << "Time,Material,Segments,Total Length,Vessels,Branch Points,Anastamoses,Active Tips,Sprouts" << endl;

    //stream3 = fopen("tracking.ang","wt");   // tracking.ang: time step, model time, total length in culture, number of branches in culture
	m_stream = fopen("out_data.ang","wt");                                        // data.ang: Store 3D coordinates of begining and end of each vessel segment
	// as well as total length of the segment}
	m_stream2 = fopen("out_vess_state.ang", "wt");
	//m_stream2 = fopen("out_vess_state.ang","wb");						// Open the stream for the vessel state data file		
	bf_stream = fopen("out_bf_state.ang","wt");						// Open the stream for the body force state data file

	time_stream = fopen("out_time.ang","wt");						// Open the stream for the time and state data file
	time_write_headers = true;										// Set the time and state data file to write the headers on its first output

	m_stream4 = fopen("out_active_tips.csv", "wt");		// active tips
	fprintf(m_stream4, "%-5s,%-12s,%-12s,%-12s\n", "State", "X", "Y", "Z");
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	//fputc(EOF, m_stream2);
	fclose(m_stream2);
	fclose(bf_stream);
	fclose(time_stream);
}

//-----------------------------------------------------------------------------
void Fileout::writeTracking(FEAngio& angio)
{
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		double l = cult->TotalVesselLength();
		SimulationTime& t = angio.CurrentSimTime();
		fprintf(m_stream3, "%-12.7f %-12.7f %-12.7f %-5i\n", t.dt, t.t, l, cult->m_num_branches);   // Write to tracking.ang
	}
}

//-----------------------------------------------------------------------------
void Fileout::closeTracking()
{
    fclose(m_stream3);                                                        // Close stream to 'tracking.ang' (stream3) 
}

//-----------------------------------------------------------------------------
void Fileout::printStatus(FEAngio& angio)
{
	//now updated to a form which can be easily consumed by excel
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		SimulationTime& t = angio.CurrentSimTime();

		
		cout << endl << "Time: " << t.t << endl;                             // Print out current time to user
		cout << "Material: " << i << endl;
		//cout << "dt: " << data.dt << endl;
		cout << "Segments     : " << cult->Segments() << endl;                             // Print out current number of segments to user
		cout << "Total Length : " << cult->TotalVesselLength() << endl;                  // Print out the current total length to user
		cout << "Vessels      : " << cult->m_num_vessel << endl;
		cout << "Branch Points: " << cult->m_num_branches << endl;                 // Print out the current number of branches to user
		cout << "Anastomoses  : " << cult->m_num_anastom << endl;            // Print out the current number of anastomoses to user
		cout << "Active tips  : " << cult->ActiveTips() << endl;
		cout << "Sprouts      : " << angio.m_pmat[i]->Sprouts() << endl;
		cout << endl;

		//Time,Material,Segments,Total Length,Vessels,Branch Points,Anastamoses,Active Tips,Sprouts
		logstream << t.t << "," << i << "," << cult->Segments() << "," << cult->TotalVesselLength() << 
			"," << cult->m_num_vessel << "," << cult->m_num_branches << "," << cult->m_num_anastom << "," 
			<< cult->ActiveTips() << "," << angio.m_pmat[i]->Sprouts() << endl;
	}
}

//-----------------------------------------------------------------------------
void Fileout::dataout(FEAngio &feangio)
{
    writeData(feangio);                                    // Create and write to 'data.ang'
}

//-----------------------------------------------------------------------------
// Will be written to data.ang
void Fileout::writeData(FEAngio &angio)
{
	fprintf(m_stream,"%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-5s %-5s\n","Time","X1","Y1","Z1","X2","Y2","Z2","Length","Vess","Label");
	
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture * cult = angio.m_pmat[i]->m_cult;
		const SegmentList& seg_list = cult->GetSegmentList();
		for (ConstSegIter it = seg_list.begin(); it != seg_list.end(); ++it)
		{
			const vec3d& r0 = it->tip(0).pos();
			const vec3d& r1 = it->tip(1).pos();
			fprintf(m_stream, "%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n", it->GetTimeOfBirth(), r0.x, r0.y, r0.z, r1.x, r1.y, r1.z, it->length(), it->m_nid, it->seed());
		}
	}
	
	fclose(m_stream);
}

//-----------------------------------------------------------------------------
void Fileout::writeNodes(FEAngio& angio)
{
 
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_nodes.ang","wt");                                       
	angio.ForEachNode([&stream2, & angio](FENode & node)
	{
		fprintf(stream2, "%-5.2i %-12.7f %-12.7f %-12.7f\n", node.GetID(), node.m_rt.x, node.m_rt.y, node.m_rt.z);
	});
                                                                   	
	fclose(stream2);                                                        
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
void Fileout::writeCollFib(FEAngio & angio, bool initial)
{
	FILE *node_stream;
	
	if (initial == true)
		node_stream = fopen("out_coll_fib_init.ang","wt");
	else
		node_stream = fopen("out_coll_fib.ang","wt");

	FEMesh * mesh = angio.GetMesh();
	for (int i = 0; i < mesh->Nodes(); i++)
	{
		FENode & node = mesh->Node(i);
		FEAngioNodeData & nd = angio.m_fe_node_data[node.GetID()];
		//TODO: remove negative oneAccessCheckAndAuditAlarm update the tests once all tests are passing
		fprintf(node_stream, "%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", node.GetID() -1, node.m_rt.x,
			node.m_rt.y, node.m_rt.z, nd.m_collfib.x, nd.m_collfib.y, nd.m_collfib.z);
	}
	fclose(node_stream);
}

//-----------------------------------------------------------------------------
void Fileout::writeECMDen(FEAngio & angio)
{
	FILE *node_stream;
	
	node_stream = fopen("out_ecm_den.ang","wt");

	FEMesh * mesh = angio.GetMesh();
	for (int i = 0; i < mesh->Nodes(); i++)
	{
		FENode & node = mesh->Node(i);
		FEAngioNodeData & nd = angio.m_fe_node_data[node.GetID()];
		//TODO: redo this once tests are passing
		fprintf(node_stream, "%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", node.GetID() -1, node.m_rt.x, node.m_rt.y,
			node.m_rt.z, nd.m_ecm_den, nd.m_ecm_den0);
	}


	fclose(node_stream);
}

//-----------------------------------------------------------------------------
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	fprintf(m_stream2,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length");  // Write column labels to out_vess_state.ang
	
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		const SegmentList& seg_list = cult->GetSegmentList();
		for (ConstSegIter it = seg_list.begin(); it != seg_list.end(); ++it)	// Iterate through all segments in frag list container (it)
		{
			const vec3d& r0 = it->tip(0).pos();
			const vec3d& r1 = it->tip(1).pos();
			fprintf(m_stream2, "%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", angio.FE_state, it->GetTimeOfBirth(), r0.x, r0.y, r0.z, r1.x, r1.y, r1.z, it->length());  // Write to out_vess_state.ang

		}
	}
}

//-----------------------------------------------------------------------------
// Save active points
void Fileout::save_active_tips(FEAngio& angio) const
{
	double t = angio.CurrentSimTime().t;
	
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture * cult = angio.m_pmat[i]->m_cult;
		const SegmentTipList& tips = cult->GetActiveSortedTipList();
		
		for (ConstTipIter it = tips.begin(); it != tips.end(); ++it)
		{
			Segment::TIP& tip = *(*it);
			if (tip.bactive)
			{
				const vec3d& r = tip.pos();
				fprintf(m_stream4, "%-5.2d,%-12.7f,%-12.7f,%-12.7f\n", angio.FE_state, r.x, r.y, r.z);
			}
		}
	}
}
void Fileout::save_timeline(FEAngio& angio)
{
#ifndef NDEBUG
	auto iter = FragmentBranching::timeline.begin();

	//consider how these shoudl be named to avoid collisions
	FILE * timeline_file = fopen("timeline.csv", "wt");
	assert(timeline_file);
	fprintf(timeline_file, "emerge_time,epoch_time,percent_of_parent,priority,callsite,branch\n");

	while (iter != FragmentBranching::timeline.end())
	{
		fprintf(timeline_file, "%-12.7f,%-12.7f,%-12.7f,%d\n", iter->emerge_time, iter->epoch_time, iter->percent_of_parent, iter->priority);
		++iter;
	}
	fclose(timeline_file);
#endif
}

//-----------------------------------------------------------------------------
// Save positions of the body forces at the current time step (This function needs to be re-written)
void Fileout::save_bdy_forces(FEAngio& angio)
{
	FEModel& fem = angio.GetFEModel();
	int NBF = fem.BodyLoads();

	fprintf(bf_stream,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length"); 

	SimulationTime& t = angio.CurrentSimTime();
	for (int i = 0; i < NBF; ++i)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
		FEParameterList& pl = pbf->GetParameterList();
		FEParam* pa = pl.Find("a");
		FEParam* prc = pl.Find("rc");

		if (pa && prc)
		{
			if (pa->value<double>() != 0.0)
				fprintf(bf_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",angio.FE_state, t.t, prc->value<vec3d>().x, prc->value<vec3d>().y, prc->value<vec3d>().z, prc->value<vec3d>().x + 1.0, prc->value<vec3d>().y + 1.0, prc->value<vec3d>().z + 1.0, 1.73); 
		}
	}
}

//-----------------------------------------------------------------------------
// Save the current time information			
void Fileout::save_time(FEAngio& angio)
{
	SimulationTime& t = angio.CurrentSimTime();
	if (time_write_headers == true){												// If this is the first time writing to out_time.ang
		fprintf(time_stream,"%-5s %-12s %-12s\n","State","Vess Time","FE Time");		// Print the column labels
		time_write_headers = false;}													// Turn off the headers flag
	
//	fprintf(time_stream,"%-5.2i %-12.7f %-12.7f\n",angio.FE_state, angio.m_time.t, ((double)angio.FE_state - 1.0)*angio.FE_time_step);	// Print out the FE state, the vessel growth model time, and the FE time
	fprintf(time_stream,"%-5.2i %-12.7f \n",angio.FE_state, t.t);	// Print out the FE state, the vessel growth model time, and the FE time
}

//-----------------------------------------------------------------------------
// Output parameter values after the simulation ends
void Fileout::output_params(FEAngio& angio)
{
	FILE *param_stream;													// Parameter output file stream                                                                                                                     
	param_stream = fopen("out_params.ang","wt");						// Output the parameter output file
	
	//fprintf(param_stream,"a = %5.5f \n", angio.m_sproutf);						// Print the sprout force magnitude
	//fprintf(param_stream,"tip range = %5.5f \n",angio.m_tip_range);				// Print the sprout force range
	//fprintf(param_stream,"phi_stiff_factor = %5.5f \n",angio.phi_stiff_factor);	// Print the displacement stiffness factor
	fprintf(param_stream,"total_body_force = %10.5i \n",angio.total_bdyf);		// Print the total number of body forces
	fclose(param_stream);												// Close the parameter output file
}
