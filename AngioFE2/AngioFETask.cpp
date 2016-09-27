// AngioFETask.cpp

#include "StdAfx.h"
#include "AngioFETask.h"
#include "FEAngio.h"
#include <iostream>
#include <FECore/FEModel.h>
using namespace std;

FEAngio* pfeangio = 0;

//-----------------------------------------------------------------------------
AngioFETask::AngioFETask(FEModel* pfem) : FECoreTask(pfem)
{
	m_pangio = 0;
}

//-----------------------------------------------------------------------------
AngioFETask::~AngioFETask(void)
{
	//TODO:find where this is deleted
	//delete m_pangio;
}

//-----------------------------------------------------------------------------
// Task initialization.
// This allocates the FEAngio object and read the angio input file. It also initializes the 
// FEAngio object as well as the FE model.
bool AngioFETask::Init(const char* inpfile)
{
	// Get the FE model
	FEModel& fem = *GetFEModel();

	cout << endl << "Angiogenesis Growth Model:" << endl << endl;

	// Create the FEAngio class
	m_pangio = new FEAngio(fem);
	pfeangio = m_pangio;

	// Read the angio3d input file
	//Filein filein;
	//if (filein.Input(inpfile, *m_pangio) == false) return false;
	
	// initialize feangio object
	if (m_pangio->Init() == false) return false;

	// all is well
	return true;
}

//-----------------------------------------------------------------------------
// This runs the actual simulation.
// This just calls the FEModel::Run function.
bool AngioFETask::Run()
{
	return GetFEModel()->Solve();
}
