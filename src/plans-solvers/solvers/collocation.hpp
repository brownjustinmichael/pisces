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
		class collocation : public solver
		{
		private:
			using plans::solvers::solver::data_in;
			using plans::solvers::solver::data_out;
			using plans::solvers::solver::element_flags;
			using plans::solvers::solver::component_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data

			mpi::messenger* messenger_ptr; //!< A pointer to the mpi messenger object

			double& timestep; //!< A double reference to the current timestep
			double *rhs_ptr; //!< A pointer to the right hand side of the equation

			int excess_0; //!< The integer number of elements to recv from edge_0
			int excess_n; //!< The integer number of elements to recv from edge_n

			double* default_matrix; //!< The double array of the non-timestep dependent matrix component

			std::vector <double> data_temp; //!< A double vector to be used in lieu of data_out for non-updating steps
			std::vector <double> matrix; //!< The left hand side of the equation in matrix form excluding the original value and the timestep mutliplication
			std::vector <double> factorized_matrix; //!< A double vector containing the factorized sum of default matrix and timestep * matrix
			std::vector <double> boundary_matrix; //!< Some memory space needed for the matrix solve
			
			std::vector <int> ns; //!< A vector containing the extents of all the elements
			std::vector <int> ipiv; //!< A vector of integers needed to calculate the factorization
			std::vector <int> bipiv; //!< A vector of integers needed to calculate the factorization
		
			std::shared_ptr <boundaries::boundary> boundary_0; //!< A shared pointer to the top boundary
			std::shared_ptr <boundaries::boundary> boundary_n; //!< A shared pointer to the bottom boundary
		
			int inner_m; //!< The number of non-overlapping grid points
			int ex_overlap_0; //!< The external overlap region on the top
			int overlap_0; //!< The total overlap region on the top
			int ex_overlap_n; //!< The external overlap region on the bottom
			int overlap_n; //!< The total overlap region on the bottom
			int lda; //!< The leading dimension of the matrices
			
		void init (mpi::messenger* i_messenger_ptr, double& i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out);

		public:
			/*!**********************************************************************
			 * @param i_messenger_ptr A pointer to the mpi messenger object
			 * \param i_timestep A double reference to the current timestep
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A reference to the data to read
			 * \param i_data_out A reference to the data to update
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
			collocation (mpi::messenger* i_messenger_ptr, double& i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out) :
			solver (i_data, i_data_out, this->get_state_in (), this->get_state ()),
			timestep (i_timestep) {
				init (i_messenger_ptr, i_timestep, i_boundary_0, i_boundary_n, i_rhs, i_data, i_data_out);
			}
			
			/*!**********************************************************************
			 * @param i_messenger_ptr A pointer to the mpi messenger object
			 * \param i_timestep A double reference to the current timestep
			 * \param i_boundary_0 A boundary factory for the top boundary
			 * \param i_boundary_n A boundary factory for the bottom boundary
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A reference to the data to read
			 * \param i_data_out A reference to the data to update
			 ************************************************************************/
			collocation (mpi::messenger* i_messenger_ptr, double& i_timestep, boundaries::boundary::factory &i_boundary_0, boundaries::boundary::factory &i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out) :
			solver (i_data, i_data_out, this->get_state_in (), this->get_state ()),
			timestep (i_timestep) {
				init (i_messenger_ptr, i_timestep, i_boundary_0.instance (i_data.get_grid (0), i_data.get_grid (1), false), i_boundary_n.instance (i_data.get_grid (0), i_data.get_grid (1), true), i_rhs, i_data, i_data_out);
			}

			virtual ~collocation () {}

			/**
			 * @copydoc solver::get_state
			 */
			int get_state () {
				return spectral_spectral;
			}
			
			/*!**********************************************************************
			 * \copydoc solver::matrix_ptr
			 ************************************************************************/
			double *matrix_ptr () {
				return &matrix [0];
			}
			
			/*!**********************************************************************
			 * \copydoc solver::factorize
			 ************************************************************************/
			void factorize ();
			
			void setup () {
				linalg::scale (m * m, 0.0, &matrix [0]);
			}
			
			/*!**********************************************************************
			 * \copydoc solver::execute
			 ************************************************************************/
			void execute ();
			
			/*!**********************************************************************
			 * \copydoc solver::factory
			 ************************************************************************/
			class factory : public plans::solvers::solver::factory
			{
			private:
				mpi::messenger *messenger_ptr; //!< A pointer to the mpi messenger object for the solver to be constructed
				double &timestep; //!< A reference to the timestep for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				boundaries::boundary::factory *boundary_factory_0; //!< A shared pointer to the top boundary for the solver to be constructed
				boundaries::boundary::factory *boundary_factory_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, double &i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n) : 
				messenger_ptr (i_messenger_ptr), 
				timestep (i_timestep), 
				boundary_0 (i_boundary_0), 
				boundary_n (i_boundary_n) {}

				/*!**********************************************************************
				 * \param i_messenger_ptr A pointer to the mpi messenger object for the solver to be constructed
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A boundary factory to the top boundary
				 * \param i_boundary_n A boundary factory to the bottom boundary
				 ************************************************************************/
				factory (mpi::messenger *i_messenger_ptr, double &i_timestep, boundaries::boundary::factory &i_boundary_0, boundaries::boundary::factory &i_boundary_n) : 
				messenger_ptr (i_messenger_ptr), 
				timestep (i_timestep), 
				boundary_factory_0 (&i_boundary_0), 
				boundary_factory_n (&i_boundary_n) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver> instance (grids::variable &i_data_in, grids::variable &i_data_out, grids::variable &i_rhs) const {
					return std::shared_ptr <plans::solvers::solver> (new collocation (messenger_ptr, timestep, boundary_factory_0 ? boundary_factory_0->instance (i_data_in.get_grids (), false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (i_data_in.get_grids (), true) : boundary_n, i_rhs.ptr (real_spectral), i_data_in, i_data_out));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: COLLOCATION_SOLVER_HPP_76DA75C5 */
