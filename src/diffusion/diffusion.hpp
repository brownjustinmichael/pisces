// 
//! \file diffusion.hpp
//  diffusion
// 
//  This is the header file for the basic diffusion classes.
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef diffusion_H_1AIIK1RA
#define diffusion_H_1AIIK1RA

#include <vector>

//! \brief Function from blas that copies a double array to another in place
//
//! \param n a pointer to an integer number of elements in x to copy to y
//! \param x the array from which the data are copied
//! \param incx a pointer to an integer spacing of elements in x
//! \param y the array to which the data are copied
//! \param incy a pointer to an integer spacing of elements in y
extern "C" void   dcopy_(int *n, double *x, int *incx, double *y, int *incy);

//! \brief Function from blas that solves a packed triangular matrix equation
//
//! \param uplo a pointer to a character, "U" for upper triangular matrix, "L" for lower triangular matrix
//! \param trans (no idea, just set this to a pointer to "N")
//! \param diag a pointer to a character, "U" for unit diagonal, "N" for anything else
//! \param n a pointer to an integer number of elements in x to copy to y
//! \param a the triangular packed array of matrix elements
//! \param x the array of right-hand-side elements in the equation which will be replaced with the resulting array
//! \param incx a pointer to an integer spacing of elements in x
extern "C" void   dtpsv_(char *uplo, char *trans, char *diag, int *n, double *a, double *x, int *incx);

namespace diffusion
{
	//! \brief The basic functional unit of the diffusion module, containing an operation and a number of data pointers
	//
	//! The plan class contains the operator and the addresses of all the relevant data arrays to operate on. Each plan need be constructed only once and can run any number of times each timestep.
	// We may decide later to remove this class from the diffusion module, as it may be generally applicable to the entire code.
	class plan
	{
	public:
		//! \brief Operate the plan operator on the data arrays specified
		//
		//! The plan class serves as a wrapper for this function. The user specifies the time step and the plan determines how to take care of the operation itself.
		//! \param timestep a double length of time over which to apply the operation
		virtual void execute (double timestep) = 0;
		// virtual ~plan ();
	};
		
	//! \brief Subclass of operation, one implementation of diffusion for data expressed as a sum of Chebyshev polynomials in 1D
	//
	//! This implementation is a simple forward time-difference scheme, solving only for the diffusion component of the equations. The class solves this implicitly with two matrix equations, one for the even and one for the odd Chebyshev polynomials. It uses the BLAS library to do so.
	class cheb_1D : public plan
	{
	private:
		double coeff; //!< a double that represents the coefficient in front of the diffusion term in the differential equation
		int n; //!< an integer number of data elements (grid points) that cheb_1D will be built to handle
		double *data_in; //!< a double pointer to the input data
		double *data_out; //!< a double pointer to the output data; if data_in == data_out, the operation is done in place
		std::vector<double> even_diffusion_matrix; //!< a 1D vector to be filled with the triangular packed matrix equation for the even Chebyshev polynomials
		std::vector<double> odd_diffusion_matrix; //!< a 1D vector to be filled with the triangular packed matrix equation for the odd Chebyshev polynomials
	public:
		//! \param i_coeff a double containing the coefficient in front of the diffusion term in the differential equation
		//! \param i_n an integer number of data elements (grid points) that cheb_1D will be built to tackle
		//! \param i_data_in a double pointer pointing to the input data
		//! \param i_data_out a double pointer pointing to the output data; if data_out == data_in or NULL, the operation is done in place
		cheb_1D (double i_coeff, int i_n, double *i_data_in, double *i_data_out = NULL);
		
		//! \brief Execute the operation on the data for a given timestep duration
		//
		//! \param timestep a double duration over which the diffusion step will happen
		void execute (double timestep);
	};	
} /* diffusion */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
