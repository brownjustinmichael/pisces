///////////////////////////////////////////////////////////////////////////////////////
//
//		Filename: advection.h
//		File type: header
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef advection_H
#define advection_H

#include <vector>

namespace advection
{

	class advec_1D
	{
	private:
	
		int n;				// number of grid points
		double *data_in;	// double pointer to input data
		double *data_out;	// double pointer to output data
		
	public:

		advec_1D (int i_n, double *i_data_in, double *i_data_out);	//constuctor initializes private members to point to input and output vectors
		virtual ~advec_1D () {}
		
		void execute (double spacestep);

	};

}

#endif