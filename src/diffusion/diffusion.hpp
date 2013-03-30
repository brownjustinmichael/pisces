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
	//! \brief Abstract base class for any diffusion solver
	// 
	//! The operation class has one virtual function, solve, which takes the timestep as a parameter and uses this to update the data.
	class operation
	{
	public:
		int n; //!< The integer number of data elements (grid points) the operator is designed to handle
		int n_data_ptrs; //!< The integer number of data arrays (e.g. input_velocity, output_temperature, etc.) the operator is designed to handle
		
		//! \param i_n an integer number of data elements (grid points) the operator should be built to handle
		//! \param i_n_data_ptrs an integer number of data arrays (e.g. input_velocity, output_temperature, etc.) the operator should be built to handle
		
		operation (int i_n, int i_n_data_ptrs) {n = i_n; n_data_ptrs = i_n_data_ptrs;}
		//! \brief A virtual function for processing the operation on the given data
		// 
		//! \param timestep a double specifying the timestep over which to update the data
		//! \param data_ptrs an array of double pointers specifying the addresses of the data arrays on which to operate
		virtual void operate (double timestep, double **data_ptrs) = 0;
		// virtual ~operation ();
	};
	
	//! \brief The basic functional unit of the diffusion module, containing an operation and a number of data pointers
	//
	//! The plan class contains the operator and the addresses of all the relevant data arrays to operate on. Each plan need be constructed only once and can run any number of times each timestep.
	// We may decide later to remove this class from the diffusion module, as it may be generally applicable to the entire code.
	class plan
	{
	private:
		operation *operation_ptr; //!< A pointer pointing to the operation for the plan
		std::vector<double *> data_ptrs; //!< A vector of double pointers that contains the addresses of all data arrays required for the execution of the plan
	public:
		//! \param i_operation_ptr a pointer to the operation specified by the plan
		//! \param i_n_data_ptrs an integer number of pointers in i_data_ptrs
		//! \param i_data_ptrs an array of double pointers to all the data arrays required for the execution of the plan
		plan (operation *i_operation_ptr, int i_n_data_ptrs, double **i_data_ptrs);
		
		//! \brief Operate the plan operator on the data arrays specified
		//
		//! The plan class serves as a wrapper for this function. The user specifies the time step and the operator determines how to take care of the operation itself.
		//! \param timestep a double length of time over which to apply the operation
		inline void execute (double timestep) {operation_ptr->operate (timestep, &data_ptrs [0]);}
		// virtual ~plan ();
	};
		
	//! \brief Subclass of operation, one implementation of diffusion for data expressed as a sum of Chebyshev polynomials in 1D
	//
	//! This implementation is a simple forward time-difference scheme, solving only for the diffusion component of the equations. The class solves this implicitly with two matrix equations, one for the even and one for the odd Chebyshev polynomials. It uses the BLAS library to do so.
	class cheb_1D : public operation
	{
	private:
		double coeff; //!< a double that represents the coefficient in front of the diffusion term in the differential equation
		std::vector<double> even_diffusion_matrix; //!< a 1D vector to be filled with the triangular packed matrix equation for the even Chebyshev polynomials
		std::vector<double> odd_diffusion_matrix; //!< a 1D vector to be filled with the triangular packed matrix equation for the odd Chebyshev polynomials
	public:
		//! \param i_n an integer number of data elements (grid points) that cheb_1D will be built to tackle
		//! \param i_coeff a double containing the coefficient in front of the diffusion term in the differential equation
		cheb_1D (int i_n, double i_coeff = 1.);
		
		//! \brief Operate on the given data for a given timestep duration
		//
		//! \param timestep a double duration over which the diffusion step will happen
		//! \param data_ptrs an array containing two pointers: the first to the original data, and the second to the location of the desired output
		void operate (double timestep, double **data_ptrs);
	};	
} /* diffusion */

#endif /* end of include guard: diffusion_H_1AIIK1RA */
