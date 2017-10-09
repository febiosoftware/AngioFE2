///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "StdAfx.h"
#include "zlib.h"

class Segment;
class FEAngio;

class Fileout
{
public:
	Fileout();
	virtual ~Fileout();
	void printStatus(FEAngio& angio);
	void save_vessel_state(FEAngio& angio);
	void save_active_tips(FEAngio& angio) const;
	void save_timeline(FEAngio& angio);
	void save_winfiber(FEAngio& angio);
	static void save_final_vessel_csv(FEAngio & angio);

private:
	std::ofstream logstream;
	FILE*	m_stream4 = 0;	// active tips
	FILE*  vessel_state_stream;
};
