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
	 * \brief \copybrief bases::solver
	 * 
	 * A LAPACK implementation of a matrix solver
	 *********************************************************************/
	class solver : public bases::solver
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
		solver (bases::element* i_element_ptr, int i_n, double& i_timestep, double& i_alpha_0, double& i_alpha_n, double *i_default_matrix, double *i_matrix, int i_name_in, int i_name_rhs, int i_name_out = null);
		
		virtual ~solver () {
			TRACE (logger, "Calling destructor.");
		}
		
		/*!*******************************************************************
		 * \copydoc bases::solver::solve ()
		 *********************************************************************/
		void execute ();
		
		virtual void update_globals (int N, double* global_matrix, double* global_rhs, int* status);
		
		virtual void update_from_globals (double* global_out);
		
		void calculate_bounds ();
		
		void send_bounds ();
		
		void recv_bounds ();
		
		void calculate_error ();
		
		void send_error ();
		
		void recv_error ();
		
		void update ();

	protected:
		double& timestep;
		double& alpha_0;
		double& alpha_n;
		double data_0, data_n;
		double prev_data_0, prev_data_n;
		
		double *rhs; //!< The double array of the right-hand-side of the matrix equation
		double* default_matrix;
		double *matrix; //!< The double matrix to be factorized
		
		std::vector <double> error;
		std::vector <double> data_temp;
		std::vector <double> factorized_matrix;
		std::vector<int> ipiv; //!< A vector of integers needed to calculate the factorization
		
		/*!*******************************************************************
		 * \copydoc bases::solver::factorize ()
		 *********************************************************************/
		void factorize ();
	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
