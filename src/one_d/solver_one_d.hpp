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

namespace bases
{
	class element;
} /* bases */

namespace one_d
{
	/*!*******************************************************************
	 * \brief A 1D implementation of a matrix solver
	 * 
	 * This matrix solver solves a matrix in each element individually and
	 * sends the results to the adjacent elements. This solver should be 
	 * iterated.
	 *********************************************************************/
	class solver : public bases::solver
	{
	public:
		/*!*******************************************************************
		 * \param i_excess_0 The integer number of excess elements on the edge_0 side
		 * \param i_excess_n The integer number of excess elements on the edge_n side
		 * \param i_timestep A double reference to the current timestep
		 * \param i_alpha_0 A double reference to the edge_0 weight
		 * \param i_alpha_n A double reference to the edge_n weight
		 * \param i_default_matrix A double array containing the 0 order collocation matrix
		 * \param i_matrix The double matrix to be factorized
		 * \param i_name_rhs The integer representation of the matrix right-hand-side
		 * \copydoc bases::solver::solver ()
		 *********************************************************************/
		solver (bases::element* i_element_ptr, int i_n, int i_excess_0, int i_excess_n, double& i_timestep, double& i_alpha_0, double& i_alpha_n, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out = null, int i_flags = 0x00);
		
		virtual ~solver () {}
		
		/*!*******************************************************************
		 * \copydoc bases::solver::execute ()
		 *********************************************************************/
		void execute ();
		
		/*!*******************************************************************
		 * \copydoc bases::solver::factorize ()
		 *********************************************************************/
		void factorize ();
		
	protected:
		double& timestep; //!< A double reference to the current timestep
		double& alpha_0; //!< A double reference to the current edge_0 weight
		double& alpha_n; //!< A double reference to the current edge_n weight

		int excess_0; //!< The integer number of elements to recv from edge_0
		int excess_n; //!< The integer number of elements to recv from edge_n
		int expected_excess_0; //!< The integer number of elements to send to edge_0
		int expected_excess_n; //!< The integer number of elements to send to edge_n
		
		double* rhs; //!< The double array of the right-hand-side of the matrix equation
		double* default_matrix; //!< The double array of the non-timestep dependent matrix component
		double* matrix; //!< The double array of the matrix component to be timestep-multiplied
		
		std::vector <double> error_0; //!< A double vector to be recved from edge_0
		std::vector <double> error_n; //!< A double vector to be recved from edge_n
		std::vector <double> out_error_0; //!< A double vector to be sent to edge_0
		std::vector <double> out_error_n; //!< A double vector to be sent to edge_n
		std::vector <double> data_temp; //!< A double vector to be used in lieu of data_out for non-updating steps
		std::vector <double> positions_0; //!< A double vector of excess positions from edge_0
		std::vector <double> positions_n; //!< A double vector of excess positions from edge_n
		std::vector <double> factorized_matrix; //!< A double vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
