#pragma once

//-----------------------------------------------------------------------------
class FEAngio;

//-----------------------------------------------------------------------------
// This class reads in the angio input file.
class Filein
{
public:
	// Constructor
	Filein();

	// Destructor
	virtual ~Filein();    

	// Read the input file
	bool Input(const char* szfile, FEAngio& angio);

private:
	// read_param - Read a parameter from the current line within the buffer
    void read_param(FEAngio &angio, char* buffer);

	// set_param - Set a certain parameter within the INPUT class based on what's read from the input file
    void set_param(FEAngio &angio, char* buffer, char* pname);		
};
