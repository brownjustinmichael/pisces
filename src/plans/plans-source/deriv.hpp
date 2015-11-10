/*!**********************************************************************
 * \file deriv.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_TWO_D_HPP_G9AN5CH6
#define SOURCE_TWO_D_HPP_G9AN5CH6

#include <omp.h>

#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"
#include "io/parameters.hpp"

/*!**********************************************************************
 * \namespace plans::source
 * 
 * \brief A namespace containing the source plan extension
 ************************************************************************/
namespace plans
{
	namespace source
	{
		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		class uniform_grad_x : public explicit_plan
		{
		private:
			using explicit_plan::coeff;
			using explicit_plan::n;
			using explicit_plan::ldn;
			using explicit_plan::m;
			using explicit_plan::dims;
			using explicit_plan::data_out;
		
			double *data_source; //!< The data pointer for the source data
			double scalar;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform_grad_x (grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : explicit_plan (i_data_in, i_data_out, i_coeff), data_source (i_data_source.ptr (real_spectral)) {
				TRACE ("Adding source...");
				scalar = acos (-1.0) * 2.0 / (i_data_in.get_grid (0) [n - 1] - i_data_in.get_grid (0) [0]);
			}
		
			virtual ~uniform_grad_x () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; i += 2)
				{
					for (int j = 0; j < m; ++j)
					{
						data_out [i * m + j] += coeff * scalar * data_source [(i + 1) * m + j];
						data_out [(i + 1) * m + j] -= coeff * scalar * data_source [i * m + j];
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public explicit_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_data_source, double i_coeff = 1.0) : explicit_plan::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new uniform_grad_x (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		class uniform_grad_z : public explicit_plan
		{
		private:
			using explicit_plan::coeff;
			using explicit_plan::n;
			using explicit_plan::ldn;
			using explicit_plan::m;
			using explicit_plan::dims;
			using explicit_plan::data_out;
		
			const double *data_source, *pos_m; //!< The data pointer for the source data
			const int ld_source;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform_grad_z (grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			explicit_plan (i_data_in, i_data_out, i_coeff), 
			data_source (i_data_source.ptr (real_spectral)),
			ld_source (i_data_source.get_ld ()) {
				TRACE ("Adding source...");
				pos_m = &(i_data_in.get_grid (1) [0]);
			}
		
			virtual ~uniform_grad_z () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; i++)
				{
					for (int j = 1; j < m - 1; ++j)
					{
						data_out [i * m + j] += coeff * (data_source [i * ld_source + j + 1] - data_source [i * ld_source + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public explicit_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_data_source, double i_coeff = 1.0) : explicit_plan::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new uniform_grad_z (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		class z_src : public real_plan
		{
		private:
			using real_plan::coeff;
			using real_plan::n;
			using real_plan::ldn;
			using real_plan::m;
			using real_plan::dims;
			using real_plan::data_out;
		
			double *data_source; //!< The data pointer for the source data
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			z_src (double *i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : real_plan (i_data_in, i_data_out, i_coeff), data_source (i_data_source) {
				TRACE ("Adding source...");
			}
		
			virtual ~z_src () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; ++i)
				{
					linalg::add_scaled (m, coeff, data_source, data_out + i * m);
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				double *data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (double *i_data_source, double i_coeff = 1.0) : real_plan::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new z_src (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

