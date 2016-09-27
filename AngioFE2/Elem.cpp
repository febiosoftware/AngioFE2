#include "StdAfx.h"
#include "Elem.h"
#include "BC.h"

FEAngioElementData::FEAngioElementData()
{
	alpha = 0.0;
	fiber_orient = vec3d();
}

FEAngioNodeData::FEAngioNodeData()
{
	m_ecm_den = 0.0;		
	m_ecm_den0 = 0.0;		

	m_collfib = vec3d();		
	m_collfib0 = vec3d();		
	vegf_conc = 0.0;		
	m_da = 1.0;
	alpha = 0.0;
	m_ntag = 0;
}