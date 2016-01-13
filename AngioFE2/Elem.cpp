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
Face::Face() : BC(false), bc_type(BC::STOP) 
{

}

//-----------------------------------------------------------------------------
BBOX Elem::GetBoundingBox()
{
	BBOX b(n1->rt);
	b.Add(n2->rt);
	b.Add(n3->rt);
	b.Add(n4->rt);
	b.Add(n5->rt);
	b.Add(n6->rt);
	b.Add(n7->rt);
	b.Add(n8->rt);
	return b;
}

//-----------------------------------------------------------------------------
void Elem::GetFace(int n, int* nf)
{
	switch (n)
	{
	case 0: nf[0] = n1->m_id; nf[1] = n2->m_id; nf[2] = n6->m_id; nf[3] = n5->m_id; break;
	case 1: nf[0] = n2->m_id; nf[1] = n4->m_id; nf[2] = n8->m_id; nf[3] = n6->m_id; break;
	case 2: nf[0] = n4->m_id; nf[1] = n3->m_id; nf[2] = n7->m_id; nf[3] = n8->m_id; break;
	case 3: nf[0] = n3->m_id; nf[1] = n1->m_id; nf[2] = n5->m_id; nf[3] = n7->m_id; break;
	case 4: nf[0] = n5->m_id; nf[1] = n6->m_id; nf[2] = n8->m_id; nf[3] = n7->m_id; break;
	case 5: nf[0] = n3->m_id; nf[1] = n4->m_id; nf[2] = n2->m_id; nf[3] = n1->m_id; break;
	}
}
