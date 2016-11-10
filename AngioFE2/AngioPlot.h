#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
class FEPlotAngioStress : public FEDomainData
{
public:
	FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioGradientCenter : public FEDomainData
{
public:
	FEPlotAngioGradientCenter(FEModel* pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioEffectiveStress : public FEDomainData
{
public:
	FEPlotAngioEffectiveStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str)  override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioCollagenFibers : public FENodeData
{
public:
	FEPlotAngioCollagenFibers(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};

class FEPlotAngioGradient : public FENodeData
{
public:
	FEPlotAngioGradient(FEModel * pfem) : FENodeData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEMesh & m, FEDataStream & a) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMDensity : public FENodeData
{
public:
	FEPlotAngioECMDensity(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMAlpha : public FENodeData
{
public:
	FEPlotAngioECMAlpha(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};