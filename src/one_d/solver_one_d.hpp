/*!***********************************************************************
 * \file one_d/solver.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-13.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOLVER_HPP_N0BUAX6H
#define SOLVER_HPP_N0BUAX6H

#include "messenger/messenger.hpp"
#include <vector>
#include <memory>
#include "plan/plan_one_d.hpp"
#include "../bases/solver.hpp"
#include "linalg/utils.hpp"

namespace one_d
{
	template <class datatype>
	class master_solver : public bases::master_solver <datatype>
	{
	private:
		int n;
		bases::grid <datatype> &grid;
		
		std::shared_ptr <bases::solver <datatype> > z_solver;
		
		std::vector <datatype> explicit_rhs_vec;
		std::vector <datatype> implicit_rhs_vec;
		std::vector <datatype> real_rhs_vec;
		
		datatype *explicit_rhs_ptr;
		datatype *implicit_rhs_ptr;
		datatype *real_rhs_ptr;
		
		using bases::master_solver <datatype>::data;
		using bases::master_solver <datatype>::element_flags;
		using bases::master_solver <datatype>::component_flags;
		
	public:
		master_solver (bases::grid <datatype> &i_grid, datatype *i_data, int *i_element_flags, int *i_component_flags) : bases::master_solver <datatype> (i_data, element_flags, component_flags), n (i_grid.get_n ()), grid (i_grid) {
			explicit_rhs_ptr = NULL;
			implicit_rhs_ptr = NULL;
			real_rhs_ptr = NULL;
		}
		
		virtual ~master_solver () {}
		
		virtual int n_dependencies () {
			if (z_solver) {
				return z_solver->n_dependencies ();
			} else {
				return 0;
			}
		}

		virtual int& get_dependency (int i) {
			return z_solver->get_dependency (i);
		}

		virtual void add_dependency (int name, int flags = 0x00) {
			if (z_solver) {
				z_solver->add_dependency (name);
			}
		}
		
		bases::grid <datatype> *grid_ptr (int index = 0) {
			return &grid;
		}
		
		datatype *rhs_ptr (int index = implicit_rhs) {
			if (index == implicit_rhs) {
				if (!implicit_rhs_ptr) {
					implicit_rhs_vec.resize (n);
					implicit_rhs_ptr = &implicit_rhs_vec [0];
				}
				return implicit_rhs_ptr;
			} else if (index == explicit_rhs) {
				if (!explicit_rhs_ptr) {
					explicit_rhs_vec.resize (n);
					explicit_rhs_ptr = &implicit_rhs_vec [0];
				}
				return explicit_rhs_ptr;
			} else if (index == real_rhs) {
				if (!real_rhs_ptr) {
					real_rhs_vec.resize (n);
					real_rhs_ptr = &implicit_rhs_vec [0];
				}
				return real_rhs_ptr;
			} else {
				return NULL;
			}
		}
		
		datatype *matrix_ptr (int index = 0) {
			return z_solver->matrix_ptr ();
		}
		
		virtual void reset () {
			if (implicit_rhs_ptr) {
				utils::scale (n, 0.0, implicit_rhs_ptr);
			}
			if (explicit_rhs_ptr) {
				utils::scale (n, 0.0, explicit_rhs_ptr);
			}
			if (real_rhs_ptr) {
				utils::scale (n, 0.0, real_rhs_ptr);
			}
		}
		
		virtual void add_solver (std::shared_ptr <bases::solver <datatype> > i_solver, int flags = 0x00) {
			z_solver = i_solver;
		}
		
		virtual std::shared_ptr <bases::solver <datatype>> get_solver (int flags = 0x00) {
			return z_solver;
		}
		
		void add_plan (const typename bases::explicit_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			bases::grid <datatype>* grids [1] = {&grid};
			bases::master_solver <datatype>::add_plan (factory.instance (grids, data, rhs_ptr (spectral_rhs), element_flags, component_flags), flags);
		}
		
		void add_plan (const typename bases::real_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			bases::grid <datatype>* grids [1] = {&grid};
			bases::master_solver <datatype>::add_plan (factory.instance (grids, data, rhs_ptr (real_rhs), element_flags, component_flags), flags);
		}
		
		void add_plan (const typename bases::implicit_plan <datatype>::factory &factory, int flags) {
			TRACE ("Adding plan...");
			bases::grid <datatype>* grids [1] = {&grid};
			datatype* matrices [1] = {matrix_ptr (0)};
			bases::master_solver <datatype>::add_plan (factory.instance (grids, matrices, data, rhs_ptr (spectral_rhs), element_flags, component_flags), flags);
		}
		
	protected:
		virtual void _factorize () {
			z_solver->factorize ();
		}
		
		virtual void _solve () {
			z_solver->execute ();
		}
	};
	
	/*!*******************************************************************
	 * \brief A 1D implementation of a matrix solver
	 * 
	 * This matrix solver solves a matrix in each element individually and
	 * sends the results to the adjacent elements. This solver should be 
	 * iterated.
	 *********************************************************************/
	template <class datatype>
	class solver : public bases::solver <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_excess_0 The integer number of excess elements on the edge_0 side
		 * \param i_excess_n The integer number of excess elements on the edge_n side
		 * \param i_timestep A datatype reference to the current timestep
		 * \param i_alpha_0 A datatype reference to the edge_0 weight
		 * \param i_alpha_n A datatype reference to the edge_n weight
		 * \param i_default_matrix A datatype array containing the 0 order collocation matrix
		 * \param i_matrix The datatype matrix to be factorized
		 * \param i_name_rhs The integer representation of the matrix right-hand-side
		 * \copydoc bases::solver <datatype>::solver ()
		 *********************************************************************/
		solver (bases::grid <datatype> &i_grid, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n, datatype *i_explicit_rhs, datatype* i_implicit_rhs, datatype *i_real_rhs, datatype* i_data, int *element_flags, int *component_flags);
		
		solver (bases::master_solver <datatype> &i_solver, utils::messenger* i_messenger_ptr, datatype& i_timestep, datatype& i_alpha_0, datatype& i_aplha_n);
	
		virtual ~solver () {
			// printf ("Destroying one_d solver\n");
		}

		datatype *matrix_ptr () {
			return &matrix [0];
		}
		
		void factorize ();
		
		void execute ();
		
		using bases::solver <datatype>::element_flags;
		using bases::solver <datatype>::component_flags;
			
	protected:
		int n, ld;
		bases::grid <datatype> &grid;
		datatype *data;
		
		utils::messenger* messenger_ptr;
		
		datatype& timestep; //!< A datatype reference to the current timestep
		datatype& alpha_0; //!< A datatype reference to the current edge_0 weight
		datatype& alpha_n; //!< A datatype reference to the current edge_n weight
		datatype value_0, value_n;
		
		const datatype* positions;
		std::vector <datatype> matrix;
		int temp_n;
		int excess_0; //!< The integer number of elements to recv from edge_0
		int excess_n; //!< The integer number of elements to recv from edge_n
		int ex_excess_0; //!< The integer number of elements to send to edge_0
		int ex_excess_n; //!< The integer number of elements to send to edge_n
		int ntop, nbot;
		datatype *implicit_rhs_vec, *explicit_rhs_vec, *real_rhs_vec;
		
		datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component
		
		std::vector <datatype> values_0, values_n;
		std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
		std::vector <datatype> positions_0; //!< A datatype vector of excess positions from edge_0
		std::vector <datatype> positions_n; //!< A datatype vector of excess positions from edge_n
		std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <datatype> boundary_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
		std::vector <datatype> previous_rhs;
		std::vector <int> ns;
		std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
		std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
	};
} /* one_d */

#endif /* end of include guard: SOLVER_HPP_N0BUAX6H */
