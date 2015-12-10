/*!**********************************************************************
 * \file plans/plans-source/uniform.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UNIFORM_HPP_G9AN5CH6
#define UNIFORM_HPP_G9AN5CH6

#include <omp.h>

#include "../explicit_plan.hpp"
#include "../real_plan.hpp"
#include "linalg/exceptions.hpp"

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
		class uniform : public explicit_plan
		{
		private:
			using explicit_plan::coeff;
			using explicit_plan::n;
			using explicit_plan::ldn;
			using explicit_plan::m;
			using explicit_plan::dims;
			using explicit_plan::data_out;
		
			double *data_source; //!< The data pointer for the source data
			bool dealias;

		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * @param i_dealias Whether to ignore the last third of the data in the horizontal
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform (grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0, bool i_dealias = false) : 
			explicit_plan (i_data_in, i_data_out, i_coeff), 
			data_source (i_data_source.ptr (real_spectral)),
			dealias (i_dealias) {
				TRACE ("Adding source...");
			}
		
			virtual ~uniform () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				linalg::matrix_add_scaled (m * dims, dealias ? ldn * 2 / 3 : ldn, coeff, data_source, data_out, m * dims, m * dims);	
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public explicit_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
				bool dealias; //!< Whether to ignore the last third of the data in the horizontal
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 * @param i_dealias Whether to ignore the last third of the data in the horizontal
				 ************************************************************************/
				factory (grids::variable &i_data_source, bool i_dealias = false, double i_coeff = 1.0) : 
				explicit_plan::factory (i_coeff), 
				data_source (i_data_source),
				dealias (i_dealias) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new uniform (data_source, i_data_in, i_data_out, coeff, dealias));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation in real-real space
		 ************************************************************************/
		class uniform_real : public real_plan
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
			 * \copydoc real_plan::real_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform_real (grids::variable &i_data_source, grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : 
			real_plan (i_data_in, i_data_out, i_coeff), 
			data_source (i_data_source.ptr (real_real)) {
				TRACE ("Adding source...");
			}
		
			virtual ~uniform_real () {}

			/*!**********************************************************************
			 * \copydoc real_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				linalg::matrix_add_scaled (m * dims, n, coeff, data_source, data_out, m * dims, m * dims);	
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public real_plan::factory
			{
			private:
				grids::variable &data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable &i_data_source, double i_coeff = 1.0) : 
				real_plan::factory (i_coeff), 
				data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new uniform_real (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a constant term to an equation
		 ************************************************************************/
		class constant : public explicit_plan
		{
		private:
			using explicit_plan::coeff;
			using explicit_plan::n;
			using explicit_plan::ldn;
			using explicit_plan::m;
			using explicit_plan::dims;
			using explicit_plan::data_out;
				
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			constant (grids::variable &i_data_in, grids::variable &i_data_out, double i_coeff = 1.0) : explicit_plan (i_data_in, i_data_out, i_coeff) {
				TRACE ("Adding constant...");
			}
		
			virtual ~constant () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; ++i)
				{
					for (int j = 0; j < m * dims; ++j)
					{
						data_out [i * m + j] += coeff;
					}
				}
			}
		
			/*!**********************************************************************
			 * \copydoc plan::factory
			 ************************************************************************/
			class factory : public explicit_plan::factory
			{
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 ************************************************************************/
				factory (double i_coeff = 1.0) : explicit_plan::factory (i_coeff) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan > _instance (double **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan > (new constant (i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: UNIFORM_HPP_G9AN5CH6 */

