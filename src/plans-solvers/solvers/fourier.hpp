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
		class fourier : public solver
		{
		private:
			using plans::solvers::solver::data_in;
			using plans::solvers::solver::data_out;
			using plans::solvers::solver::element_flags;
			using plans::solvers::solver::component_flags;
			
			int n; //!< The horizontal extent of the data
			int ldn; //!< The horizontal extent of the data array
			int m; //!< The vertical extent of the data
			
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
			
			const double *pos_m; //!< A pointer to the array of positions in the vertical

			void init (double& i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out);
			
		public:
			/*!**********************************************************************
			 * \param i_timestep A double reference to the current timestep
			 * \param i_boundary_0 A shared pointer to the top boundary object
			 * \param i_boundary_n A shared pointer to the bottom boundary object
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A reference to the data to read
			 * \param i_data_out A reference to the data to update
			 ************************************************************************/
			fourier (double& i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out) :
			solver (i_data, i_data_out, this->get_state_in (), this->get_state ()),
			timestep (i_timestep) {
				init (i_timestep, i_boundary_0, i_boundary_n, i_rhs, i_data, i_data_out);
			}
			
			/*!**********************************************************************
			 * \param i_timestep A double reference to the current timestep
			 * \param i_boundary_0 A boundary factory for the top boundary
			 * \param i_boundary_n A boundary factory for the bottom boundary
			 * \param i_rhs A pointer to the right hand side of the equation
			 * \param i_data A reference to the data to read
			 * \param i_data_out A reference to the data to update
			 ************************************************************************/
			fourier (double& i_timestep, typename boundaries::boundary::factory &i_boundary_0, typename boundaries::boundary::factory &i_boundary_n, double *i_rhs, grids::variable &i_data, grids::variable &i_data_out) :
			solver (i_data, i_data_out, this->get_state_in (), this->get_state ()),
			timestep (i_timestep) {
				init (i_timestep, i_boundary_0.instance (i_data.get_grid (0), i_data.get_grid (1), false), i_boundary_n.instance (i_data.get_grid (0), i_data.get_grid (1), true), i_rhs, i_data, i_data_out);
			}
			
			virtual ~fourier () {}
			
			/**
			 * @copydoc solver::get_state
			 */
			int get_state () {
				return real_spectral;
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
			
			/**
			 * @copydoc plan::setup
			 */
			void setup () {
				linalg::scale (m * ldn, 0.0, &matrix [0]);
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
				double &timestep; //!< The timestep reference for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_0; //!< A shared pointer to the top boundary for the solver to be constructed
				std::shared_ptr <boundaries::boundary> boundary_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				typename boundaries::boundary::factory *boundary_factory_0; //!< A shared pointer to the top boundary for the solver to be constructed
				typename boundaries::boundary::factory *boundary_factory_n; //!< A shared pointer to the bottom boundary for the solver to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A shared pointer to the top boundary for the solver to be constructed
				 * \param i_boundary_n A shared pointer to the bottom boundary for the solver to be constructed
				 ************************************************************************/
				factory (double &i_timestep, std::shared_ptr <boundaries::boundary> i_boundary_0, std::shared_ptr <boundaries::boundary> i_boundary_n) : 
				timestep (i_timestep), 
				boundary_0 (i_boundary_0), 
				boundary_n (i_boundary_n) {}

				/*!**********************************************************************
				 * \param i_timestep A reference to the timestep for the solver to be constructed
				 * \param i_boundary_0 A reference to a boundary factory to generate the top boundary
				 * \param i_boundary_n A reference to a boundary factory to generate the bottom boundary
				 ************************************************************************/
				factory (double &i_timestep, typename boundaries::boundary::factory &i_boundary_0, typename boundaries::boundary::factory &i_boundary_n) : 
				timestep (i_timestep), 
				boundary_factory_0 (&i_boundary_0), 
				boundary_factory_n (&i_boundary_n) {}

				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc solver::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::solvers::solver> instance (grids::variable &i_data_in, grids::variable &i_data_out, grids::variable &i_rhs) const {
					return std::shared_ptr <plans::solvers::solver> (new fourier (timestep, boundary_factory_0 ? boundary_factory_0->instance (i_data_in.get_grids (), false) : boundary_0, boundary_factory_n ? boundary_factory_n->instance (i_data_in.get_grids (), true) : boundary_n, i_rhs.ptr (real_spectral), i_data_in, i_data_out));
				}
			};
		};
	} /* solvers */
} /* plans */

#endif /* end of include guard: FOURIER_SOLVER_HPP_6DE77A57 */
