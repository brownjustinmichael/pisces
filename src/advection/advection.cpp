///////////////////////////////////////////////////////////////////////////////////////
//
//		Filename: advection.cpp
//		File type: source
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "advection.h"

namespace advection
{
	
	advec_1D::advec_1D (int i_n, double *i_data_in, double *i_data_out)
	{
		n = i_n;				// initialize number of grid points
		data_in = i_data_in;	// initialize pointer to input data
		data_out = i_data_out;	// initialize pointer to output data
	}

	void advec_1D::execute(double spacestep)
	{
		int i;
		double inv_dx = 1/spacestep;

		*data_out = inv_dx * ( *(data_in + 1) - *data_in );										// make initial forward difference step
		for (i = 1; i < n - 1; i++)
		{
				*(data_out + i) = .5 * inv_dx * ( *(data_in + i + 1) - *(data_in + i - 1) );	// centered differencing scheme
		}
		*data_out = inv_dx * ( *(data_in + n - 1) - *(data_in + n - 2) );						// make backward difference step for final grid point
	}

}