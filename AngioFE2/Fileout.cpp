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
#include <unordered_set>


using namespace std;

//-----------------------------------------------------------------------------
Fileout::Fileout()
{
    logstream.open("out_log.csv");
	//write the headers
	logstream << "Time,Material,Segments,Total Length,Vessels,Branch Points,Anastamoses,Active Tips,Sprouts" << endl;

	vessel_state_stream = fopen("out_vess_state.ang2" , "wb");//check the parameters consider setting the compression level
	unsigned int magic = 0xfdb97531;
	unsigned int version = 0;
	fwrite(&magic, sizeof(unsigned int), 1, vessel_state_stream);
	fwrite(&version, sizeof(unsigned int), 1, vessel_state_stream);


	m_stream4 = fopen("out_active_tips.csv", "wt");		// active tips
	fprintf(m_stream4, "%-5s,%-12s,%-12s,%-12s\n", "State", "X", "Y", "Z");
}

//-----------------------------------------------------------------------------
Fileout::~Fileout()
{
    logstream.close();
	fclose(vessel_state_stream);
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
// Save microvessel position at the current time point
void Fileout::save_vessel_state(FEAngio& angio)
{
	unsigned int segcount = 0;
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		auto seg_list = cult->PerStepSegments();
		segcount += seg_list.size();
	}
	//write segcount and time
	fwrite(&segcount, sizeof(unsigned int), 1, vessel_state_stream);
	SimulationTime st = angio.CurrentSimTime();
	float rtime = static_cast<float>(st.t);
	fwrite(&rtime, sizeof(float), 1, vessel_state_stream);

	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		auto seg_list = cult->PerStepSegments();
		for (auto it = seg_list.begin(); it != seg_list.end(); ++it)	// Iterate through all segments in frag list container (it)
		{
			vec3d r0 = angio.ReferenceCoords((*it)->tip(0).pt);
			vec3d r1 = angio.ReferenceCoords((*it)->tip(1).pt);
			//consider checking fwrites return values
			float cur = static_cast<float>(r0.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r0.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);

			cur = static_cast<float>(r1.x);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.y);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
			cur = static_cast<float>(r1.z);
			fwrite(&cur, sizeof(float), 1, vessel_state_stream);
		}
	}

	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture* cult = angio.m_pmat[i]->m_cult;
		cult->ClearPerStepSegments();
	}

}

//-----------------------------------------------------------------------------
// Save active points
void Fileout::save_active_tips(FEAngio& angio) const
{	
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

void Fileout::save_final_vessel_csv(FEAngio & angio)
{
	auto iter = FragmentBranching::timeline.begin();

	//consider how these shoudl be named to avoid collisions
	FILE * final_vessel_file = fopen("final_vessels.csv", "wt");
	assert(final_vessel_file);
	fprintf(final_vessel_file, "x0,y0,z0,x1,y1,z1,start time\n");

	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture * cult = angio.m_pmat[i]->m_cult;
		auto segs = cult->GetSegmentList();

		for (auto it = segs.begin(); it != segs.end(); ++it)
		{
			auto seg = (*it);
			fprintf(final_vessel_file, "%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f,%-12.7f\n", seg.tip(0).pt.r.x, seg.tip(0).pt.r.y, seg.tip(0).pt.r.z, seg.tip(1).pt.r.x, seg.tip(1).pt.r.y, seg.tip(1).pt.r.z, seg.GetTimeOfBirth());
		}
	}

	fclose(final_vessel_file);
}

void Fileout::save_winfiber(FEAngio& angio)
{
	FILE * winfiber_file = fopen("final_state.mv3d", "wt");
	assert(winfiber_file);
	fprintf(winfiber_file, "#a file for WInFiber3d Generated by AngioFE\n");
	
	//maintain a way to get consistent pointers
	std::map<int, std::vector<Segment*>> vessels;
	int num_segs = 0;
	for (size_t i = 0; i < angio.m_pmat.size(); i++)
	{
		Culture * cult = angio.m_pmat[i]->m_cult;
		num_segs += cult->m_frag.size();
		for(auto j = cult->m_frag.begin();j != cult->m_frag.end();++j)
		{
			if(vessels.count(j->m_nvessel))
			{
				vessels[j->m_nvessel].push_back(&(*j));
			}
			else
			{
				std::vector<Segment*> sl;
				vessels[j->m_nvessel] = sl;
				vessels[j->m_nvessel].push_back(&(*j));
			}
		}
	}

	int winfib_counter = 0;

	for(auto i = vessels.begin(); i != vessels.end(); ++i)
	{
		std::map<int, bool> seg_processed;
		//pick a segment and visit all connected segments
		for(auto j = 0; j < i->second.size(); ++j)
		{
			if(!seg_processed.count(i->second.at(j)->m_nid))
			{
				Segment * cur = i->second.at(j);
				std::unordered_set<int> visited;
				visited.insert(cur->m_nid);
				while(cur->tip(0).connected && visited.count(cur->tip(0).connected->m_nid) == 0)
				{
					cur = cur->tip(0).connected;
					visited.insert(cur->m_nid);
				}
				fprintf(winfiber_file, "%d\t%f\t%f\t%f\t%f\n", winfib_counter, cur->tip_c(0).pt.r.x, cur->tip_c(0).pt.r.y, cur->tip_c(0).pt.r.z, cur->length());
				
				while(cur->tip(1).connected)
				{
					//write cur to file
					fprintf(winfiber_file, "%d\t%f\t%f\t%f\t%f\n", winfib_counter, cur->tip_c(1).pt.r.x, cur->tip_c(1).pt.r.y, cur->tip_c(1).pt.r.z, cur->length());
					seg_processed[cur->m_nid] = true;
					cur = cur->tip(1).connected;
				}
				seg_processed[cur->m_nid] = true;
				fprintf(winfiber_file, "%d\t%f\t%f\t%f\t%f\n\n", winfib_counter, cur->tip_c(1).pt.r.x, cur->tip_c(1).pt.r.y, cur->tip_c(1).pt.r.z, cur->length());
				winfib_counter++;
				//assert(seg_processed.count(i->second.at(j)->m_nid));
			}
		}
	}

	fclose(winfiber_file);
}


