///////////////////////////////////////////////////////////////////////////////////////
//
//		! \file one_d/advection.cpp
//		File type: source
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include "advection_one_d.h"

namespace one_d
{
	
	advec::advec (int i_n, double *i_tmstp_ptr, double i_c, double *i_data_in, double *i_data_out)
	{
		int i;

		n = i_n;				// initialize number of grid points
		c = i_c;
		tmstp_ptr = i_tmstp_ptr;
		fac = (n - 1) / acos(-1.0);
		data_in = i_data_in;	// initialize pointer to input data
		data_out = i_data_out;	// initialize pointer to output data

		sin_vals.resize(n, 0.0);
		for (i = 1; i < n-1; i++)
		{
			sin_vals [i] = fac / sin (i/fac);
			// sin_vals [i] = 0;
			// sin_vals [i] = fac;
			// sin_vals [i] = fac / (cos ((n - i - 2)/fac) - cos ((n - i)/fac));
		}
	}
	
	advec::advec (int i_n, double *i_tmstp_ptr, double i_c, double& i_data_in, double& i_data_out)
	{
		int i;

		n = i_n;				// initialize number of grid points
		c = i_c;
		tmstp_ptr = i_tmstp_ptr;
		fac = (n - 1) / acos(-1.0);
		data_in = &i_data_in;	// initialize pointer to input data
		data_out = &i_data_out;	// initialize pointer to output data

		sin_vals.resize(n, 0.0);
		sin_vals [0] = fac / sin (0.5 / fac);
		sin_vals [n - 1] = fac / sin ((n - 1.5) / fac);
		for (i = 1; i < n-1; i++)
		{
			sin_vals [i] = fac / sin (i/fac);
			// sin_vals [i] = 0;
			// sin_vals [i] = fac;
			// sin_vals [i] = fac / (cos ((n - i - 2)/fac) - cos ((n - i)/fac));
		}
	}

	void advec::execute()
	{
		int i;

		*data_out += (*tmstp_ptr) * c * sin_vals [0] * ( *(data_in + 1) - *data_in );																						// Impose left boundary condition
		for (i = 1; i < n - 1; i++)
		{
				*(data_out + i) += 0.5 * (*tmstp_ptr) * c * sin_vals [i] * ( *(data_in + i + 1) - *(data_in + i - 1) );			// centered differencing scheme
		}
		*(data_out + n - 1) += (*tmstp_ptr) * c * sin_vals [n - 1] * ( *(data_in + n - 1) - *(data_in + n - 2));																			// Impose right boundary condition
	}

}