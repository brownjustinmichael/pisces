/*!**********************************************************************
 * \file variable_diffusion.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VARIABLE_DIFFUSION_HPP_E1A9847D
#define VARIABLE_DIFFUSION_HPP_E1A9847D

#include "linalg/utils.hpp"

#include "../real_plan.hpp"
#include "../implicit_plan.hpp"

namespace plans
{
	namespace diffusion
	{
		/**
		 * @brief An implicit plan designed to include diffusion with an additional multiplier
		 * @details The precise form of the term is grad dot ($source) grad ($data). This usually occurs in the inclusion of density to the equations, though that is not the only case where terms like this appear. The execution of this plan does use a bit of a hack for stability: the right hand side generated is that from a full timestep ago.
		 */
		
		class variable_diffusion : public implicit_plan
		{
		protected:
			using implicit_plan::n;
			using implicit_plan::ldn;
			using implicit_plan::m;
			using implicit_plan::grid_n;
			using implicit_plan::grid_m;
			using implicit_plan::data_in;
			using implicit_plan::data_out;
			using implicit_plan::coeff;
			using implicit_plan::matrix_n;
			using implicit_plan::matrix_m;
			
			using implicit_plan::element_flags;
			using implicit_plan::component_flags;

			double *data_source; //!< A pointer to the source data
			double *new_matrix; //!< A pointer to the diffusion matrix
			double *current; //!< A pointer to the current values of the diffusion
			double *bg_deriv; //!< A pointer to the derivative of the background source
			double *bg_deriv2; //!<  pointer to the second derivative of the background source
			std::vector<double> bg_deriv_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<double> bg_deriv2_vec; //!< A vector containing the derivative of the background diffusion
			std::vector<double> new_matrix_vec; //!< A vector containing the implicit diffusion matrix
			std::vector<double> bg_val; //!< A vector containing the average background source
			std::vector<double> current_vec; //!< A vector containing the current data
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary

		public:
			/**
			 * @param matrix_n A pointer to the horizontal matrix
			 * @param matrix_m A pointer to the vertical matrix
			 * @param i_data_source A reference to the source value that enters inside the divergence
			 * @param i_data_in A reference to the input data
			 * @param i_data_out A reference to the output data
			 * @param i_coeff The coefficient by which to multiply the results
			 */
			variable_diffusion (double *matrix_n, double *matrix_m, grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0);
			
			virtual ~variable_diffusion () {}

			/**
			 * @copydoc plan::setup
			 */
			void setup ();

			/**
			 * @copydoc plan::execute
			 */
			virtual void execute ();

			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public implicit_plan::factory
			{
			private:
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 ************************************************************************/
				factory (grids::variable &i_data_source, double i_coeff = 1.0) : implicit_plan::factory (i_coeff), data_source (i_data_source) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const;
			};
		};

		/*!**********************************************************************
		 * \brief A plan to enact a diffusion coefficient linearly dependent on a variable
		 ************************************************************************/
		
		class linear : public real_plan
		{
		protected:
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::grid_n;
			using real_plan::grid_m;
			using real_plan::data_in;
			using real_plan::data_out;
			
			using real_plan::element_flags;
			using real_plan::component_flags;
			
			int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
			int count; //!< The current count of evaluations
			
			double coeff; //!< The base coefficient to multiply the source
			double min; //!< The minimum value of the coefficient
			double *data_source; //!< A pointer to the source data
			double pioL2; //!< The value pi / L^2 for speed
			
			std::vector <double> x1_vec; //!< A vector for an intermediary step
			std::vector <double> x2_vec; //!< A vector for an intermediary step
			std::vector <double> z1_vec; //!< A vector for an intermediary step
			std::vector <double> z2_vec; //!< A vector for an intermediary step
			std::vector <double> bg_state_vec; //!< A vector to contain the background state for the diffusion coefficient
			std::vector <double> bg_deriv_vec; //!< A vector ot contain the derivative of the background state for the diffusion coefficient

			double *x1_ptr; //!< A pointer to the x1 vector
			double *x2_ptr; //!< A pointer to the x2 vector
			double *z1_ptr; //!< A pointer to the z1 vector
			double *z2_ptr; //!< A pointer to the x2 vector
			double *bg_state; //!< A pointer to the bg_state vector
			double *bg_diff; //!< A pointer to the background diffusion coefficient
			double *bg_deriv; //!< A pointer to the bg_deriv vector
			const double *pos_m; //!< A pointer to the vertical position, for speed
			const double *pos_n; //!< A pointer to the horizontal positions, for speed
			double *oodx; //!< An array of one over the differential of x, centered in the cell
			double *oodx2; //!< An array of one over the differential of x, centered at the boundary
			double *oodz; //!< An array of one over the differential of z, centered in the cell
			double *oodz2; //!< An array of one over the differential of z, centered at the boundary
			
		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			linear (double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0);
			
			virtual ~linear () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			void execute ();

			void reset_source (double *new_source);
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			protected:
				using real_plan::factory::coeff;
				double min; //!< The minimum value for the coefficient
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				double *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (double i_coeff, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every) : real_plan::factory (i_coeff), min (i_min), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const;
			};
		};

		class tanh : public linear
		{
		protected:
			using linear::n;
			using linear::ldn;
			using linear::m;
			using linear::grid_n;
			using linear::grid_m;
			using linear::data_in;
			using linear::data_out;
			
			using linear::element_flags;
			using linear::component_flags;
			
			double source_length;
			double source_zero;

			std::vector <double> src_vec;
			double *src_ptr;
			double *data_orig; //!< A pointer to the source data

		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			tanh (double i_source_length, double i_source_zero, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0);
			
			virtual ~tanh () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			void execute ();
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			protected:
				using real_plan::factory::coeff;
				double min; //!< The minimum value for the coefficient
				double source_length;
				double source_zero;
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				double *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (double i_coeff, double i_source_length, double i_source_zero, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every) : real_plan::factory (i_coeff), min (i_min), source_length (i_source_length), source_zero (i_source_zero), data_source (i_data_source), bg_diff (i_bg_diff), bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const;
			};
		};

		class ra_tanh : public linear
		{
		protected:
			using linear::n;
			using linear::ldn;
			using linear::m;
			using linear::grid_n;
			using linear::grid_m;
			using linear::data_in;
			using linear::data_out;
			
			using linear::element_flags;
			using linear::component_flags;
			
			double source_length;
			double source_zero;
			double a, b, c, d, e;

			std::vector <double> src_vec;
			double *src_ptr;
			double *data_orig; //!< A pointer to the source data

		public:
			/*!**********************************************************************
			 * \copydoc real_plan::real_plan
			 * 
			 * @param i_min The minimum value of the coefficient
			 * \param i_data_source A pointer to the source data
			 * @param i_bg_diff The background diffusion coefficient
			 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
			 ************************************************************************/
			ra_tanh (double i_source_length, double i_source_zero, double chi, double stiffness, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0);
			
			virtual ~ra_tanh () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			void execute ();
			
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			protected:
				using real_plan::factory::coeff;
				double min; //!< The minimum value for the coefficient
				double source_length;
				double source_zero;
				double phi;
				double stiffness;
				grids::variable &data_source; //!< The source data pointer for the plan to be constructed
				double *bg_diff; //!< The background diffusion coefficient
				int bg_every; //!< The number of evaluations between which to recalculate the background diffusion
				
			public:
				/*!**********************************************************************
				 * \param i_coeff The base coefficient for the plan to be constructed
				 * \param i_data_source The source data pointer for the plan to be constructed
				 * @param i_min The minimum value of the data for the plan to be constructed
				 * @param i_bg_diff The background diffusion coefficient
				 * @param i_bg_every The number of evaluations between which to recalculate the background diffusion
				 ************************************************************************/
				factory (double i_coeff, double i_source_length, double i_source_zero, double i_phi, double i_stiffness, double i_min, grids::variable &i_data_source, double *i_bg_diff, int i_bg_every) :
					real_plan::factory (i_coeff),
					min (i_min),
					source_length (i_source_length),
					source_zero (i_source_zero),
					phi(i_phi),
					stiffness(i_stiffness),
					data_source (i_data_source),
					bg_diff (i_bg_diff),
					bg_every (i_bg_every) {}
				
				virtual ~factory () {}
				
				/*!**********************************************************************
				 * \copydoc plan::factory::instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const;
			};
		};
	} /* diffusion */
} /* plans */

#endif /* end of include guard: VARIABLE_DIFFUSION_HPP_E1A9847D */
