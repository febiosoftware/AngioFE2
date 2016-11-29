///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <fstream>
#include <vector>

class Segment;
class FEAngio;

using namespace std;

class Fileout
{
public:
	Fileout();
	virtual ~Fileout();
	void writeTracking(FEAngio& angio);
	void closeTracking();
	void printStatus(FEAngio& angio);
	void dataout(FEAngio &feangio);
    void writeData(FEAngio &feangio);
	void writeNodes(FEAngio& angio);
	void writeEconn(FEAngio& angio);
	void writeCollFib(FEAngio & angio, bool initial);
	void writeECMDen(FEAngio & angio);
	void save_vessel_state(FEAngio& angio);
	void save_active_tips(FEAngio& angio);
	void save_bdy_forces(FEAngio& angio);
	void save_time(FEAngio& angio);
	void output_params(FEAngio& angio);

private:
	ofstream logstream;
	FILE*	m_stream;	// Open stream to 'data.ang' (stream)
	FILE*	m_stream2;	// vessel state 
    FILE*	m_stream3;
	FILE*	bf_stream;
	FILE*	m_stream4;	// active tips

	FILE *time_stream;
	bool time_write_headers;
};
