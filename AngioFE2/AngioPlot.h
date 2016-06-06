#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
class FEPlotAngioStress : public FEDomainData
{
public:
	FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str);
};

//-----------------------------------------------------------------------------
class FEPlotAngioEffectiveStress : public FEDomainData
{
public:
	FEPlotAngioEffectiveStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str);
};

//-----------------------------------------------------------------------------
class FEPlotAngioCollagenFibers : public FENodeData
{
public:
	FEPlotAngioCollagenFibers(FEModel* pfem) : FENodeData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMDensity : public FENodeData
{
public:
	FEPlotAngioECMDensity(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMAlpha : public FENodeData
{
public:
	FEPlotAngioECMAlpha(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a);
};