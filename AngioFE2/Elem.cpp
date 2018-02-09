#include "StdAfx.h"
#include "Elem.h"
#include "BC.h"

FEAngioElementData::FEAngioElementData()
{
	alpha = 0.0;
	flags = 0;
}

FEAngioNodeData::FEAngioNodeData()
{
	alpha = 0.0;
	m_ntag = 0;
}