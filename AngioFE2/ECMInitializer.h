#pragma once
#include "FiberManager.h"

class FEAngioMaterialBase;

//needs to set the intitial ecm density and anisotropy
//if the density has not yet been modified it will be 0.0
class ECMInitializer
{
public:
	virtual ~ECMInitializer() {}
	virtual void seedECMDensity(FEAngioMaterialBase* mat) = 0;
	virtual bool overwrite() { return true; }
	virtual void updateECMdensity(FEAngioMaterialBase* mat);
};
class ECMInitializerConstant : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterialBase* mat) override;
};
class ECMInitializerSpecified : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterialBase* mat) override;
};
class ECMInitializerNoOverwrite : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterialBase* mat) override;
	bool overwrite() override { return false; }
};
