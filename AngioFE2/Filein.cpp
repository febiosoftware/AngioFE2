///////////////////////////////////////////////////////////////////////
// Filein.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Filein.h"
#include <iostream>
#include <fstream>
#include "FEAngio.h"

using namespace std;

//-----------------------------------------------------------------------------
Filein::Filein()
{
    
}

//-----------------------------------------------------------------------------
Filein::~Filein()
{

}

//-----------------------------------------------------------------------------
bool Filein::Input(const char* szfile, FEAngio& angio)
{
	FILE* fp = fopen(szfile, "rt");
	if (fp == 0) return false;

    const int buff_max = 100;										// Maximum size of the buffer
    char buffer[buff_max];											// Create the buffer
    char first_char;												// First character of the current line in the buffer
    char str_type = 'n';											// Set default string type to 'no line'
    
    while (!feof(fp))												// Read until you reach the end of the file...
    {
        str_type = 'n';												// Reset the default string type to 'no line'
		fgets(buffer, 255, fp);
                
        first_char = buffer[0];										// Determine the first character of the line in the buffer
       
        if (first_char == 37)										// If first character is a '%' symbol...
             str_type = 'c';										// ... then set string type to 'Comment'
             else if (first_char == 62)								// If first character is a '>' symbol...
             str_type = 'p';										// ... then set string type to 'Param'
                      
        if (str_type == 'p')
            read_param(angio, buffer);								// If the line describes a parameter, then read in that parameter
    }

	fclose(fp);

	return true;
}

//-----------------------------------------------------------------------------
void Filein::read_param(FEAngio& angio, char* buffer)
{
    char pname[20] = {0};											// Create string for the parameter name
    
    sscanf(buffer,"%*s %s %*f",&pname);								// Scan the buffer for the name of the parameter
    
    set_param(angio,buffer,pname);									// Set the parameter based on the name
    
    return;
}

//-----------------------------------------------------------------------------
void Filein::set_param(FEAngio& angio, char* buffer, char* pname)
{
    //// Parameters for angio3d (In paratheneses, string identifying the parameter and an example value):
	// Branching Probability (brnch_ch 0.1)
	if (!strcmp(pname,"brnch_ch")){
		sscanf(buffer,"%*s %*s %f",&angio.data.branch_chance);
        return;}

    //  Matrix conditions (matx_cnd 0) random
	if (!strcmp(pname,"matx_cnd")){
        sscanf(buffer,"%*s %*s %i",&angio.grid.load_cond);
        return;}
    
	// Initial matrix density (matx_den 3.0) mg/mL
    if (!strcmp(pname,"matx_den")){
        sscanf(buffer,"%*s %*s %lf",&angio.grid.coll_den);
        return;}  
        
	// Number of initial fragments (nfrag 70, based on 30K frags/mL)
    if (!strcmp(pname,"nfrag")){
        sscanf(buffer,"%*s %*s %i",&angio.data.NFRAGS);
        return;}        
    
	// End of culture period (max_time 6.0) days
    if (!strcmp(pname,"max_time")){
        sscanf(buffer,"%*s %*s %lf",&angio.data.maxt);
        return;}  
    
	// Initial time step (dt 0.25) days
    if (!strcmp(pname,"dt")){
        sscanf(buffer,"%*s %*s %lf",&angio.data.dt);
        return;}  
    
	// Anastomosis distance (anst_dst 25.0) um 
    if (!strcmp(pname,"anst_dst")){
		sscanf(buffer,"%*s %*s %lf",&angio.data.anast_dist);
        return;}  
        
    // Segment length adjustment scale (lngth_adj 1.0) 
	if (!strcmp(pname,"lngth_adj")){
		sscanf(buffer,"%*s %*s %lf",&angio.data.length_adjust);
        return;}
    
	// Number of nodes in x-direction for autogrid (xnodes 7)
    if (!strcmp(pname,"xnodes")){
        sscanf(buffer,"%*s %*s %i",&angio.grid.xnodes);
        return;}  
    
	// Number of nodes in y-direction for autogrid (ynodes 7)
    if (!strcmp(pname,"ynodes")){
        sscanf(buffer,"%*s %*s %i",&angio.grid.ynodes);
        return;} 
    
	// Number of nodes in z-direction for autogrid (znodes 3)
    if (!strcmp(pname,"znodes")){
        sscanf(buffer,"%*s %*s %i",&angio.grid.znodes);
        return;} 
    
	// Alternative way of specifying the number of nodes in each direction.  Reads three int in succession indicating the number of nodes in the x-, y-, and z-
	// direction, respectively (num_nodes 7 7 3)
    if (!strcmp(pname,"num_nodes")){
        sscanf(buffer,"%*s %*s %i %i %i",&angio.grid.xnodes,&angio.grid.ynodes,&angio.grid.znodes);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the x-direction (xrange 0 300) um
    if (!strcmp(pname,"xrange")){
		sscanf(buffer,"%*s %*s %lf %lf",&angio.grid.xrange[0],&angio.grid.xrange[1]);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the y-direction (yrange 0 300) um
    if (!strcmp(pname,"yrange")){
		sscanf(buffer,"%*s %*s %lf %lf",&angio.grid.yrange[0],&angio.grid.yrange[1]);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the z-direction (zrange 0 300) um
    if (!strcmp(pname,"zrange")){
		sscanf(buffer,"%*s %*s %lf %lf",&angio.grid.zrange[0],&angio.grid.zrange[1]);
        return;}
    
	// Specify if a composite material model is being using 
	if (!strcmp(pname,"comp_mat")){
        sscanf(buffer,"%*s %*s %i",&angio.comp_mat);
        return;}
    
	// Specify if a composite material model is being using 
	if (!strcmp(pname,"sproutf")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_sproutf);
        return;}

	// Specify if a composite material model is being using 
	if (!strcmp(pname,"tip_range")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_tip_range);
        return;}

	// Specify if a composite material model is being using 
	if (!strcmp(pname,"sub_cyc")){
		sscanf(buffer,"%*s %*s %i",&angio.m_sub_cycles);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"gelbc")){
        sscanf(buffer,"%*s %*s %c",&angio.m_cgelbc);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"stiff_fact")){
		sscanf(buffer,"%*s %*s %lf",&angio.phi_stiff_factor);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"sprout_sphere")){
		sscanf(buffer,"%*s %*s %i",&angio.m_bsp_sphere);
        return;}

		// Specify which type of boundary condition at the front edge of the gel
	if (!strcmp(pname,"front_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.frontbc);
        return;}
	
	// Specify which type of boundary condition at the right edge of the gel
	if (!strcmp(pname,"right_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.rightbc);
        return;}
	
	// Specify which type of boundary condition at the back edge of the gel
	if (!strcmp(pname,"back_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.backbc);
        return;}
	
	// Specify which type of boundary condition at the left edge of the gel
	if (!strcmp(pname,"left_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.leftbc);
        return;}
	
	// Specify which type of boundary condition at the bottom edge of the gel
	if (!strcmp(pname,"bottom_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.bottombc);
        return;}
	
	// Specify which type of boundary condition at the top edge of the gel
	if (!strcmp(pname,"top_bc")){
        sscanf(buffer,"%*s %*s %c",&angio.grid.topbc);
        return;}

	// Specify the location of the x symmetry panel
	if (!strcmp(pname,"Sx")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_Sx);
        return;}

	// Specify the location of the y symmetry panel
	if (!strcmp(pname,"Sy")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_Sy);
        return;}

	// Specify the location of the z symmetry panel
	if (!strcmp(pname,"Sz")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_Sz);
        return;}

	// Specify the seed number for the random generator
	if (!strcmp(pname,"rseed")){
        sscanf(buffer,"%*s %*s %i",&angio.m_irseed);
        return;}

	// Specify if branching is on or off
	if (!strcmp(pname,"branch")){
		int n;
        sscanf(buffer,"%*s %*s %i",&n); angio.yes_branching = (n != 0);
        return;}
	
	// Specify if anastomosis is on or off
	if (!strcmp(pname,"anast")){
		int n;
        sscanf(buffer,"%*s %*s %i",&n); angio.yes_anast = (n != 0);
        return;}

	// Specify the factor for the directional sprout force
	if (!strcmp(pname,"spfact")){
        sscanf(buffer,"%*s %*s %lf",&angio.m_spfactor);
        return;}

	// Specify to 'flatten' fibers in z
	if (!strcmp(pname,"zfibflat")){
        sscanf(buffer,"%*s %*s %lf",&angio.grid.m_bzfibflat);
        return;}

	// Read in weights for determing the direction of growth
    	if (!strcmp(pname,"gweights")){
        sscanf(buffer,"%*s %*s %lf %lf",&angio.cult.W[0],&angio.cult.W[3]);
        return;}

    return;
}
