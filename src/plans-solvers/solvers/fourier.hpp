/*!**********************************************************************
 * \file fourier.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FOURIER_SOLVER_HPP_6DE77A57
#define FOURIER_SOLVER_HPP_6DE77A57

#include "mpi/messenger.hpp"

#include "../solver.hpp"
#include "../boundaries/boundary.hpp"

namespace plans
{
	namespace solvers
	{
		/*!**********************************************************************
		 * \brief Solve the data in the horizontal direction by fourier solve
		 * 
		 * This class is a simple implementation of the diagonal solve needed to solve the fourier components in the horizontal. If doing the full solve proves to be less complex than expected, this class will not be needed.
		 ************************************************************************/
		template <class datatype>
		class fourier : public solver <datatype>
		{
		private:
			using plans::solvers::solver <datatype>::element_flags;
			using plans::solvers::solver <datatype>::component_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data
			datatype *data; //!< A pointer to the data
			
			datatype& timestep; //!< A datatype reference to the current timestep
			datatype *rhs_ptr; //!< A pointer to the right hand side of the equation
			
			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n
			
			datatype* default_matrix; //!< The datatype array of the non-timestep dependent matrix component
			
			std::vector <datatype> data_temp; //!< A datatype vector to be used in lieu of data_out for non-updating steps
			std::vector <datatype> matrix; //!< The left hand side of the equation in matrix form excluding the original value and the timestep mutliplication
			std::vector <datatype> factorized_matrix; //!< A datatype vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <datatype> boundary_matrix; //!< Some memory space needed for the matrix solve
			
			std::vector <int> ns; //!< A vector containing the extents of all the elements
			std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
		
			std::shared_ptr <boundaries::boundary <datatype>> boundary_0; //!< A shared pointer to the top boundary
			std::shared_ptr <boundaries::boundary <datatype>> boundary_n; //!< A shared pointer to the bottom boundary
		
			int inner_m; //!< The number of non-overlapping grid points
			int ex_overlap_0; //!< The external overlap region on the top
			int overlap_0; //!< The total overlap region on the top
			int ex_overlap_n; //!< The external overlap region on the bottom
			int overlap_n; //!< The total overlap region on the bottom
			int lda; //!< The leading dimension of the matrices
			
			const datatype *pos_m; //!< A pointer to the array of positions in the vertical
			
		public:
			/*!**********************************************************************
			 * \copydoc solver::solver
			 * 
			 * \param i_grid_n The horizontal grid object
			 * \param i_grid_m The vertical grid object
			 * \param timestep A datatype reference to the current timestep
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A pointer to the data
			 ************************************************************************/
			fourier (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, datatype& i_timestep, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags);
			
			virtual ~fourier () {}
			
			/*!**********************************************************************
			 * \copydoc solver::matrix_ptr
			 ************************************************************************/
			datatype *matrix_ptr () {
				return &matrix [0];
			}
			
			/*!**********************************************************************
			 * \copydoc solver::factorize
			 ************************************************************************/
			void factorize ();
			
			/*!**********************************************************************
			 * \copydoc solver::execute
			 ************************************************************************/
			void execute ();
			
			/*!**********************************************************************
			 * \copydoc solver::factory
			 ************************************************************************/
			class factory : public plans::solvers::solver <datatype>::factory
			{
			private:
				datatype &timestep; //!< The timestep reference for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 ************************************************************************/
				factory (datatype &i_timestep, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n) : timestep (i_timestep), boundary_0 (i_boundary_0), boundary_n (i_boundary_n) {}
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver <datatype>> instance (grids::grid <datatype> **grids, datatype *i_data, datatype *i_rhs, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					return std::shared_ptr <plans::solvers::solver <datatype>> (new fourier (*grids [0], *grids [1], timestep, boundary_0, boundary_n, i_rhs, i_data, i_element_flags, i_component_flags));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: FOURIER_SOLVER_HPP_6DE77A57 */
