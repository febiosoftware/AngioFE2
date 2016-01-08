#include "stdafx.h"
#include "Elem.h"

//-----------------------------------------------------------------------------
Node::Node() : theta(0.), eta(0.), theta0(0.), eta0(0.), ecm_den(0.), ecm_den0(0.), id(0), updated(false) 
{
}

//-----------------------------------------------------------------------------
Node::Node(const Node& n)
{
	r0 = n.r0;
	rt = n.rt;

	theta = n.theta;
	theta0 = n.theta0;

	eta = n.eta;
	eta0 = n.eta0;

	ecm_den = n.ecm_den;
	ecm_den0 = n.ecm_den0;

	id = n.id;

	updated = n.updated;

	ecm_den_grad = n.ecm_den_grad;
	u = n.u;

	collfib = n.collfib;
	collfib0 = n.collfib0;

	ecm_den_store = n.ecm_den_store;
	ecm_fibril_store = n.ecm_fibril_store;
}

//-----------------------------------------------------------------------------
void Node::operator = (const Node& n)
{
	r0 = n.r0;
	rt = n.rt;

	theta = n.theta;
	theta0 = n.theta0;

	eta = n.eta;
	eta0 = n.eta0;

	ecm_den = n.ecm_den;
	ecm_den0 = n.ecm_den0;

	id = n.id;

	updated = n.updated;

	ecm_den_grad = n.ecm_den_grad;
	u = n.u;

	collfib = n.collfib;
	collfib0 = n.collfib0;

	ecm_den_store = n.ecm_den_store;
	ecm_fibril_store = n.ecm_fibril_store;
}

//-----------------------------------------------------------------------------
double Elem::bb_xmin()
{
    double exmin = n1->rt.x;
    if ((*n2).rt.x < exmin) exmin = (*n2).rt.x;
    if ((*n3).rt.x < exmin) exmin = (*n3).rt.x;
    if ((*n4).rt.x < exmin) exmin = (*n4).rt.x;
    if ((*n5).rt.x < exmin) exmin = (*n5).rt.x;
    if ((*n6).rt.x < exmin) exmin = (*n6).rt.x;
    if ((*n7).rt.x < exmin) exmin = (*n7).rt.x;
	if ((*n8).rt.x < exmin) exmin = (*n8).rt.x;

    return exmin;
}

double Elem::bb_xmax()
{
    double exmax = n1->rt.x;
    if ((*n2).rt.x > exmax) exmax = (*n2).rt.x;
    if ((*n3).rt.x > exmax) exmax = (*n3).rt.x;
    if ((*n4).rt.x > exmax) exmax = (*n4).rt.x;
    if ((*n5).rt.x > exmax) exmax = (*n5).rt.x;
    if ((*n6).rt.x > exmax) exmax = (*n6).rt.x;
    if ((*n7).rt.x > exmax) exmax = (*n7).rt.x;
	if ((*n8).rt.x > exmax) exmax = (*n8).rt.x;

	return exmax;
}

double Elem::bb_ymin()
{
	double eymin = n1->rt.y;
    if ((*n2).rt.y < eymin) eymin = (*n2).rt.y;
    if ((*n3).rt.y < eymin) eymin = (*n3).rt.y;
    if ((*n4).rt.y < eymin) eymin = (*n4).rt.y;
    if ((*n5).rt.y < eymin) eymin = (*n5).rt.y;
    if ((*n6).rt.y < eymin) eymin = (*n6).rt.y;
    if ((*n7).rt.y < eymin) eymin = (*n7).rt.y;
	if ((*n8).rt.y < eymin) eymin = (*n8).rt.y;

    return eymin;
}

double Elem::bb_ymax()
{
    double eymax = n1->rt.y;
    if ((*n2).rt.y > eymax) eymax = (*n2).rt.y;
    if ((*n3).rt.y > eymax) eymax = (*n3).rt.y;
    if ((*n4).rt.y > eymax) eymax = (*n4).rt.y;
    if ((*n5).rt.y > eymax) eymax = (*n5).rt.y;
    if ((*n6).rt.y > eymax) eymax = (*n6).rt.y;
    if ((*n7).rt.y > eymax) eymax = (*n7).rt.y;
	if ((*n8).rt.y > eymax) eymax = (*n8).rt.y;

	return eymax;
}
 
double Elem::bb_zmin()
{
 	double ezmin = n1->rt.z;
    if ((*n2).rt.z < ezmin) ezmin = (*n2).rt.z;
    if ((*n3).rt.z < ezmin) ezmin = (*n3).rt.z;
    if ((*n4).rt.z < ezmin) ezmin = (*n4).rt.z;
    if ((*n5).rt.z < ezmin) ezmin = (*n5).rt.z;
    if ((*n6).rt.z < ezmin) ezmin = (*n6).rt.z;
    if ((*n7).rt.z < ezmin) ezmin = (*n7).rt.z;
	if ((*n8).rt.z < ezmin) ezmin = (*n8).rt.z;

    return ezmin;
}

double Elem::bb_zmax()
{
    double ezmax = n1->rt.z;
    if ((*n2).rt.z > ezmax) ezmax = (*n2).rt.z;
    if ((*n3).rt.z > ezmax) ezmax = (*n3).rt.z;
    if ((*n4).rt.z > ezmax) ezmax = (*n4).rt.z;
    if ((*n5).rt.z > ezmax) ezmax = (*n5).rt.z;
    if ((*n6).rt.z > ezmax) ezmax = (*n6).rt.z;
    if ((*n7).rt.z > ezmax) ezmax = (*n7).rt.z;
	if ((*n8).rt.z > ezmax) ezmax = (*n8).rt.z;

	return ezmax;
}       
