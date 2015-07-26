/*!**********************************************************************
 * \file collocation.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_SOLVER_HPP_76DA75C5
#define COLLOCATION_SOLVER_HPP_76DA75C5

#include "linalg-block/block_solver.hpp"

#include "../boundaries/boundary.hpp"
#include "../solver.hpp"

namespace plans
{
	namespace solvers
	{
		/*!**********************************************************************
		 * \brief A solver class that solves by collocation in the vertical
		 * 
		 * This is the main solver of the code. It used both the messenger and boundary classes.
		 ************************************************************************/
		template <class datatype>
		class collocation : public solver <datatype>
		{
		private:
			using plans::solvers::solver <datatype>::data_in;
			using plans::solvers::solver <datatype>::data_out;
			using plans::solvers::solver <datatype>::element_flags;
			using plans::solvers::solver <datatype>::component_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data

			mpi::messenger* messenger_ptr; //!< A pointer to the mpi messenger object

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
			
		public:
			/*!**********************************************************************
			 * \copydoc solver::solver
			 * 
			 * \param i_grid_n The horizontal grid object
			 * \param i_grid_m The vertical grid object
			 * \param i_messenger_ptr A pointer to the mpi messenger object
			 * \param i_timestep A datatype reference to the current timestep
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A pointer to the data
			 * 
			 * The collocation matrix is set up as 
			 * 
			 * 0 0 boundary row for above element       0 0
			 * 0 0 interpolating row for above element  0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [matrix                            ] 0 0
			 * 0 0 [boundary row for this element     ] 0 0
			 * 0 0 [interpolating row for this element] 0 0
			 * 0 0 interpolating row for below element  0 0
			 * 0 0 boundary row for below element       0 0
			 ************************************************************************/
			collocation (mpi::messenger* i_messenger_ptr, datatype& i_timestep, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n, datatype *i_rhs, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out);
			
			collocation (mpi::messenger* i_messenger_ptr, datatype& i_timestep, typename boundaries::boundary <datatype>::factory &i_boundary_0, typename boundaries::boundary <datatype>::factory &i_boundary_n, datatype *i_rhs, grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out) : collocation (i_messenger_ptr, i_timestep, i_boundary_0.instance (i_data.get_grid (0), i_data.get_grid (1), false), i_boundary_n.instance (i_data.get_grid (0), i_data.get_grid (1), true), i_rhs, i_data, i_data_out) {}

			virtual ~collocation () {}

			int get_state () {
				return spectral_spectral;
			}
			
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
			
			void setup () {
				DEBUG ("FIXING");
				linalg::scale (m * m, 0.0, &matrix [0]);
			}
			
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
				mpi::messenger *messenger_ptr; //!< A pointer to the mpi messenger object for the solver to be constructed
				datatype &timestep; //!< A reference to the timestep for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary <datatype>> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				typename boundaries::boundary <datatype>::factory *boundary_factory_0; //!< A shared pointer to the top boundary for the solver to be constructed
				typename boundaries::boundary <datatype>::factory *boundary_factory_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, datatype &i_timestep, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_0, std::shared_ptr <boundaries::boundary <datatype>> i_boundary_n) : messenger_ptr (i_messenger_ptr), timestep (i_timestep), boundary_0 (i_boundary_0), boundary_n (i_boundary_n) {}

				factory (mpi::messenger *i_messenger_ptr, datatype &i_timestep, typename boundaries::boundary <datatype>::factory &i_boundary_0, typename boundaries::boundary <datatype>::factory &i_boundary_n) : messenger_ptr (i_messenger_ptr), timestep (i_timestep), boundary_factory_0 (&i_boundary_0), boundary_factory_n (&i_boundary_n) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver <datatype>> instance (grids::variable <datatype> &i_data, grids::variable <datatype> &i_data_out, grids::variable <datatype> &i_rhs) const {
					return std::shared_ptr <plans::solvers::solver <datatype>> (new collocation (messenger_ptr, timestep, boundary_factory_0 ? boundary_factory_0->instance (i_data.get_grids (), false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (i_data.get_grids (), true) : boundary_n, i_rhs.ptr (real_spectral), i_data, i_data_out));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: COLLOCATION_SOLVER_HPP_76DA75C5 */
