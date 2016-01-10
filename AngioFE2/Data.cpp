#include "stdafx.h"
#include "Grid.h"
#include "Data.h"
#include "angio3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Data::Data()
{
	// Initial time (Default: 0 days)
    t = 0.0;

    dt_a = 0.0637;
    dt_b = 9.0957;
    dt_c = 2.6073;
    n = 1.0;

    a1 = -1.2653;                                               // a1, a2, a3 - Parameters for branching curve
	a2 = 1.535;
	a3 = 0.1108;

    a = 1900.0;                                              // a, b, x0, y0 - Parameters for growth curve
    b = 1.4549;
    x0 = 4.9474;
    y0 = -19.1278;

    m_d = y0 + a/(1+pow(E,x0/b));                                 // d - Initial value of growth curve (t = 0)

	vessel_width = 7;                                           // vessel_width - Diameter of microvessels (Default: 7 um)
	num_anastom = 0;                                            // num_anastom - Initialize anastomoses counter
	num_branches = 0;                                           // num_branches - Initialize branching counter
	num_zdead = 0;
	
	branch = false;                                             // branch - Set branching indicator to 'false'

	nsegs = 0;                                                  // nsegs - Initialize segment counter
    num_vessel = NFRAGS - 1;                                    // num_vessel - Initialize vessel counter
        
	vess_length = m_d;
	old_length = m_d;

	num_lines = -1;

	NFRAGS = 3;
    maxt = 0.0;
	dt = 0.25;
	length_adjust = 1.0;
	anast_dist = 75.0;
	branch_chance = 0.1;
}


Data::~Data()                                                   // Destructor for DATA object
{

}

void Data::init_data(Grid &grid)
{
	dx = (grid.xrange[1] - grid.xrange[0])/(grid.xnodes-1);   // dx - Spatial step size in x direction, range of x divided by the number of nodes
	dy = (grid.yrange[1] - grid.yrange[0])/(grid.ynodes-1);   // dy - Spatial step size in y direction, range of y divided by the number of nodes
	dz = (grid.zrange[1] - grid.zrange[0])/(grid.znodes-1);   // dz - Spatial step size in z direction, range of z divided by the number of nodes
}
