#pragma once

class FEAngioMaterial;

//needs to set the intitial ecm density and anisotropy
//if the density has not yet been modified it will be 0.0
class ECMInitializer
{
public:
	virtual ~ECMInitializer() {}
	virtual void seedECMDensity(FEAngioMaterial * mat) = 0;
	virtual bool overwrite() { return true; }
	virtual void updateECMdensity(FEAngioMaterial * mat);
};
class ECMInitializerConstant : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
};
class ECMInitializerSpecified : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
};
class ECMInitializerNoOverwrite : public ECMInitializer
{
	void seedECMDensity(FEAngioMaterial * mat) override;
	bool overwrite() override { return false; }
};
