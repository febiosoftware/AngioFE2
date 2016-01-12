///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <list>
#include <time.h>
#include <fstream>
#include <vector>

class Grid;
class Segment;
class FEAngio;

using namespace std;

class Fileout
{
public:
	Fileout();
	virtual ~Fileout();
	
public:
	void writeTracking(FEAngio& angio);
	void closeTracking();
	void printStatus(FEAngio& angio);
	void dataout(FEAngio &feangio);
    void writeData(FEAngio &feangio);
	void writeNodes(FEAngio& angio);
	void writeEconn(FEAngio& angio);
	void writeCollFib(Grid &grid, bool initial);
	void writeECMDen(Grid &grid);
	void writeECMDenGrad(Grid &grid);
	void writeBC(Grid &grid);
	void printtime(FEAngio& angio);
	void printrandseed(int randseed);
	void writeECMDenStore(Grid &grid);
	void writeECMFibrilStore(Grid &grid);
	void save_vessel_state(FEAngio& angio);
	void save_bdy_forces(FEAngio& angio);
	void save_time(FEAngio& angio);

private:
	ofstream logstream;
    FILE*	m_stream3;      
	FILE*	m_stream;     // Open stream to 'data.ang' (stream)
	FILE*	m_stream2; 
	FILE*	bf_stream;

	FILE *time_stream;
	bool time_write_headers;
};
