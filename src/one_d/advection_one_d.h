///////////////////////////////////////////////////////////////////////////////////////
//
//		! \file one_d/advection.h
//		File type: header
//		Author: Ryan Moll (ryan.d.moll@gmail.com)
//		Date created: April 4, 2013
//		Description: Contains methods for explicit handling of spatial derivatives
//
///////////////////////////////////////////////////////////////////////////////////////

#ifndef advection_H
#define advection_H

#include <memory>
#include <vector>
#include "../bases/plan.hpp"

namespace one_d
{

	class advec : public bases::plan
	{
	private:
	
		int n;							// number of grid points
		double c;
		double fac;
		double &timestep;
		double *data_in;				// double pointer to input data
		double *data_out;				// double pointer to output data
		std::vector<double> sin_vals;	// double vector of sine values

	public:

		advec (int i_n, double& i_timestep, double i_c, double *i_data_in, double *i_data_out);	//constuctor initializes private members to point to input and output vectors
		advec (int i_n, double & i_timestep, double i_c, double& i_data_in, double& i_data_out);	//constuctor initializes private members to point to input and output vectors
		
		virtual ~advec () {}
		
		void execute ();
	};

}

#endif