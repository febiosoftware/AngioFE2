#include "stdafx.h"
#include "Elem.h"
#include "BC.h"

//-----------------------------------------------------------------------------
Node::Node()
{
	m_id = 0;
	m_ecm_den = 0.0;
	m_ecm_den0 = 0.0;
	m_ntag = 0;
}

//-----------------------------------------------------------------------------
Face::Face()
{
	bc_type = BC::STOP;
	m_nelem = -1;
}

//-----------------------------------------------------------------------------
Elem::Elem()
{
	elem_num = -1;
	
	m_volume = 0.0;
	m_volume0 = 0.0;
	alpha = 0.0;

	m_pface[0] = 0;
	m_pface[1] = 0;
	m_pface[2] = 0;
	m_pface[3] = 0;
	m_pface[4] = 0;
	m_pface[5] = 0;
}

//-----------------------------------------------------------------------------
BBOX Elem::GetBoundingBox()
{
	BBOX b(m_pnode[0]->rt);
	b.Add(m_pnode[1]->rt);
	b.Add(m_pnode[2]->rt);
	b.Add(m_pnode[3]->rt);
	b.Add(m_pnode[4]->rt);
	b.Add(m_pnode[5]->rt);
	b.Add(m_pnode[6]->rt);
	b.Add(m_pnode[7]->rt);
	return b;
}

//-----------------------------------------------------------------------------
void Elem::GetFace(int n, int* nf)
{
	switch (n)
	{
	case 0: nf[0] = m_pnode[0]->m_id; nf[1] = m_pnode[1]->m_id; nf[2] = m_pnode[5]->m_id; nf[3] = m_pnode[4]->m_id; break;
	case 1: nf[0] = m_pnode[1]->m_id; nf[1] = m_pnode[3]->m_id; nf[2] = m_pnode[7]->m_id; nf[3] = m_pnode[5]->m_id; break;
	case 2: nf[0] = m_pnode[3]->m_id; nf[1] = m_pnode[2]->m_id; nf[2] = m_pnode[6]->m_id; nf[3] = m_pnode[7]->m_id; break;
	case 3: nf[0] = m_pnode[2]->m_id; nf[1] = m_pnode[0]->m_id; nf[2] = m_pnode[4]->m_id; nf[3] = m_pnode[6]->m_id; break;
	case 4: nf[0] = m_pnode[4]->m_id; nf[1] = m_pnode[5]->m_id; nf[2] = m_pnode[7]->m_id; nf[3] = m_pnode[6]->m_id; break;
	case 5: nf[0] = m_pnode[2]->m_id; nf[1] = m_pnode[3]->m_id; nf[2] = m_pnode[1]->m_id; nf[3] = m_pnode[0]->m_id; break;
	}
}
