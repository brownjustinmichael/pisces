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
#include "../utils/utils.hpp"

namespace one_d
{
	
	advec::advec (int i_n, double& i_timestep, double i_c, double *i_data_in, double *i_data_out, std::shared_ptr<bases::collocation_grid> i_grid) : timestep (i_timestep)
	{
		int i;
		grid = i_grid;
		n = i_n;				// initialize number of grid points
		c = i_c;
		timestep = i_timestep;
		fac = (n - 1) / acos(-1.0);
		data_in = i_data_in;	// initialize pointer to input data
		data_out = i_data_out;	// initialize pointer to output data

		sin_vals.resize(n, 0.0);

		sin_vals [0] = fac / sin (.5/fac);
		for (i = 1; i < (n-1); i++)
		{
			sin_vals [i] = fac / sin (i/fac);
			// sin_vals [i] = 0;
			// sin_vals [i] = fac;
			// sin_vals [i] = fac / (cos ((n - i - 2)/fac) - cos ((n - i)/fac));
		}
		sin_vals [n-1] = fac / sin ((n-1.5)/fac);
	}
	
	advec::advec (int i_n, double& i_timestep, double i_c, double &i_data_in, double &i_data_out, std::shared_ptr<bases::collocation_grid> i_grid) : timestep (i_timestep)
	{
		int i;
		grid = i_grid;
		n = i_n;				// initialize number of grid points
		c = i_c;
		timestep = i_timestep;
		fac = (n - 1) / acos(-1.0);
		data_in = &i_data_in;	// initialize pointer to input data
		data_out = &i_data_out;	// initialize pointer to output data

		sin_vals.resize(n, 0.0);

		sin_vals [0] = fac / sin (.5/fac);
		for (i = 1; i < (n-1); i++)
		{
			sin_vals [i] = fac / sin (i/fac);
			// sin_vals [i] = 0;
			// sin_vals [i] = fac;
			//sin_vals [i] = fac / (cos ((n - i - 2)/fac) - cos ((n - i)/fac));
		}
		sin_vals [n-1] = fac / sin ((n-1.5)/fac);
	}

	void advec::execute()
	{
		double scalar = c * timestep;

		// utils::matrix_vector_multiply (n - 2, n, scalar, grid->get_data (1) + 1, data_in, 1.0, data_out + 1, n);
		
		
		*data_out += (timestep) * c * sin_vals [0] * ( *(data_in + 1) - *data_in);																											// Impose left boundary condition
		for (int i = 1; i < (n-1); i++)
		{
				*(data_out + i) += 0.5 * (timestep) * c * sin_vals [i] * ( *(data_in + i + 1) - *(data_in + i - 1) );			// centered differencing scheme
		}
		*(data_out + n - 1) += (timestep) * c * sin_vals [n-1] * ( *(data_in + n - 1) - *(data_in + n - 2));																								// Impose right boundary condition
		
	}

}