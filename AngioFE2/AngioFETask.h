#pragma once
#include <FECore/FECoreKernel.h>
#include <FECore/FECoreTask.h>

//-----------------------------------------------------------------------------
class FEAngio;

//-----------------------------------------------------------------------------
class AngioFETask : public FECoreTask
{
public:
	explicit AngioFETask(FEModel* pfem);
	~AngioFETask(void);

	bool Init(const char* szfile) override;

	bool Run() override;

private:
	FEAngio*	m_pangio=nullptr;
};
