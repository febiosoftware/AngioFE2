#pragma once
#include <FECore/FECoreKernel.h>
#include <FECore/FECoreTask.h>

//-----------------------------------------------------------------------------
class FEAngio;

//-----------------------------------------------------------------------------
class AngioFETask : public FECoreTask
{
public:
	AngioFETask(FEModel* pfem);
	~AngioFETask(void);

	bool Init(const char* szfile);

	bool Run();

private:
	FEAngio*	m_pangio;
};
