#pragma once
#include <FECore/FEPlotData.h>

//-----------------------------------------------------------------------------
class FEPlotAngioStress : public FEDomainData
{
public:
	explicit FEPlotAngioStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixStress : public FEDomainData
{
public:
	explicit FEPlotMatrixStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixViscoStress : public FEDomainData
{
public:
	explicit FEPlotMatrixViscoStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotMatrixElasticStress : public FEDomainData
{
public:
	explicit FEPlotMatrixElasticStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotMatrixElastic_m_Q : public FEDomainData
{
public:
	explicit FEPlotMatrixElastic_m_Q(FEModel* pfem) : FEDomainData(PLT_MAT3F, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotVesselStress : public FEDomainData
{
public:
	explicit FEPlotVesselStress(FEModel* pfem) : FEDomainData(PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotVesselWeight : public FEDomainData
{
public:
	explicit FEPlotVesselWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixWeight : public FEDomainData
{
public:
	explicit FEPlotMatrixWeight(FEModel* pfem) : FEDomainData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};
//-----------------------------------------------------------------------------
class FEPlotMatrixTangent : public FEDomainData
{
public:
	explicit FEPlotMatrixTangent(FEModel* pfem) : FEDomainData(PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioGradientCenter : public FEDomainData
{
public:
	explicit FEPlotAngioGradientCenter(FEModel* pfem);
	bool Save(FEDomain& d, FEDataStream& str) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioMaterialHop : public FEDomainData
{
public:
	explicit FEPlotAngioMaterialHop(FEModel* pfem);
	bool Save(FEDomain& d, FEDataStream& str) override;
};


//-----------------------------------------------------------------------------
class FEPlotAngioSegmentBadGrowth : public FEDomainData
{
public:
	explicit FEPlotAngioSegmentBadGrowth(FEModel* pfem);
	bool Save(FEDomain& d, FEDataStream& str) override;
};

class FEPlotAngioGradient : public FENodeData
{
public:
	explicit FEPlotAngioGradient(FEModel * pfem) : FENodeData(PLT_VEC3F, FMT_ITEM){}
	bool Save(FEMesh & m, FEDataStream & a) override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMDensity : public FENodeData
{
public:
	explicit FEPlotAngioECMDensity(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};

//-----------------------------------------------------------------------------
class FEPlotAngioECMAlpha : public FENodeData
{
public:
	explicit FEPlotAngioECMAlpha(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM){}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};

/*
//-----------------------------------------------------------------------------
class FEPlotAngioNodeMQ : public FENodeData
{
public:
	FEPlotAngioNodeMQ(FEModel* pfem) : FENodeData(PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEMesh& m, FEDataStream& a)  override;
};
*/

//stopgap plot variables 
class FEPlotMatrixConectrationGradient : public FEDomainData
{
public:
	explicit FEPlotMatrixConectrationGradient(FEModel * pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

class FEPlotMatrixSBMConectration : public FEDomainData
{
public:
	explicit FEPlotMatrixSBMConectration(FEModel * pfem) : FEDomainData(PLT_VEC3F, FMT_ITEM) {}
	bool Save(FEDomain& d, FEDataStream& str) override;
};

	