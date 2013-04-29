/*!***********************************************************************
 * \file one_d/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_N0BUAX6H
#define SOLVER_HPP_N0BUAX6H

#include <vector>
#include <memory>
#include "../bases/solver.hpp"
#include "../utils/utils.hpp"

namespace one_d
{
	/*!*******************************************************************
	 * \brief \copybrief bases::solver
	 * 
	 * A LAPACK implementation of a matrix solver
	 *********************************************************************/
	class lapack_solver : public bases::solver
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data_in The double array of input
		 * \param i_rhs The double array of the right-hand-side of the matrix equation
		 * \param i_matrix The double matrix to be factorized
		 * \param i_data_out The double array of output
		 * \copydoc bases::solver::solver ()
		 *********************************************************************/
		lapack_solver (int i_n, double *i_data_in, double *i_rhs, double *i_matrix, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : bases::solver (i_flags_ptr, i_logger) {
			init (i_n, i_data_in, i_rhs, i_matrix, i_data_out, i_flags_ptr);
		};
		
		/*
			TODO change location of matrix in arguments
		*/
		
		lapack_solver (int i_n, double& i_data_in, double& i_rhs, double *i_matrix, double& i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) : bases::solver (i_flags_ptr, i_logger) {
			init (i_n, &i_data_in, &i_rhs, i_matrix, &i_data_out, i_flags_ptr);
		};
		
		lapack_solver (int i_n, double& i_data_in, double& i_rhs, double *i_matrix, int *i_flags_ptr = NULL, int i_logger = -1) : bases::solver (i_flags_ptr, i_logger) {
			init (i_n, &i_data_in, &i_rhs, i_matrix);
		};
		
		virtual ~lapack_solver () {}
		
		/*!*******************************************************************
		 * \copydoc bases::solver::factorize ()
		 *********************************************************************/
		void factorize ();
		
		/*!*******************************************************************
		 * \copydoc bases::solver::solve ()
		 *********************************************************************/
		void solve ();

	private:
		int n; //!< The integer number of elements in the data
		
		double *data_in; //!< The double array of input
		double *rhs; //!< The double array of the right-hand-side of the matrix equation
		double *matrix; //!< The double matrix to be factorized
		double *data_out; //!< The double array of output
		
		std::vector<int> ipiv; //!< A vector of integers needed to calculate the factorization
		
		void init (int i_n, double *i_data_in, double *i_rhs, double *i_matrix, double *i_data_out = NULL, int *i_flags_ptr = NULL) {
			n = i_n;
			data_in = i_data_in;
			rhs = i_rhs;
			matrix = i_matrix;
			if (i_data_out) {
				data_out = i_data_out;
			} else {
				data_out = i_data_in;
			}
			ipiv.resize (i_n, 0);
		}
	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
