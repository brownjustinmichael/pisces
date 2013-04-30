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
#include "../bases/collocation.hpp"
#include "../bases/plan.hpp"

extern "C" void dgemv_ (char *trans, int *m, int *n, double *alpha, double *a, int *lda, double *x, int *incx, double *beta, double *y, int *incy);

namespace one_d
{

	class advec : public bases::plan
	{
	private:
	
		int n;							// number of grid points
		double c;
		double fac;
		double& timestep;
		double *data_in;				// double pointer to input data
		double *data_out;				// double pointer to output data
		std::vector<double> sin_vals;	// double vector of sine values
		std::shared_ptr<bases::collocation_grid> grid;

	public:

		advec (int i_n, double& i_timestep, double i_c, double *i_data_in, double *i_data_out, std::shared_ptr<bases::collocation_grid> i_grid);	//constuctor initializes private members to point to input and output vectors
		advec (int i_n, double& i_timestep, double i_c, double &i_data_in, double &i_data_out, std::shared_ptr<bases::collocation_grid> i_grid);	//constuctor initializes private members to point to input and output vectors
		
		virtual ~advec () {}
		
		void execute ();
	};

}

#endif