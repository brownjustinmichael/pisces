/*!**********************************************************************
 * \file plans/plans-source/source.hpp
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
		template <class datatype>
		class uniform : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_out;
		
			datatype *data_source; //!< The data pointer for the source data
			bool dealias;

		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform (grids::variable <datatype> &i_data_source, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0, bool i_dealias = false) : 
			explicit_plan <datatype> (i_data_in, i_data_out, i_coeff), 
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &data_source; //!< The data source to be used when constructing the plan
				bool dealias;
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_source, bool i_dealias = false, datatype i_coeff = 1.0) : 
				explicit_plan <datatype>::factory (i_coeff), 
				data_source (i_data_source),
				dealias (i_dealias) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform <datatype> (data_source, i_data_in, i_data_out, coeff, dealias));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class uniform_grad_x : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_out;
		
			datatype *data_source; //!< The data pointer for the source data
			datatype scalar;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			uniform_grad_x (grids::variable <datatype> &i_data_source, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : explicit_plan <datatype> (i_data_in, i_data_out, i_coeff), data_source (i_data_source.ptr (real_spectral)) {
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_source, datatype i_coeff = 1.0) : explicit_plan <datatype>::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform_grad_x <datatype> (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class uniform_grad_z : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_out;
		
			const datatype *data_source, *pos_m; //!< The data pointer for the source data
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
			uniform_grad_z (grids::variable <datatype> &i_data_source, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
			explicit_plan <datatype> (i_data_in, i_data_out, i_coeff), 
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				grids::variable <datatype> &data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_source, datatype i_coeff = 1.0) : explicit_plan <datatype>::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new uniform_grad_z <datatype> (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class z_src : public real_plan <datatype>
		{
		private:
			using real_plan <datatype>::coeff;
			using real_plan <datatype>::n;
			using real_plan <datatype>::ldn;
			using real_plan <datatype>::m;
			using real_plan <datatype>::dims;
			using real_plan <datatype>::data_out;
		
			datatype *data_source; //!< The data pointer for the source data
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			z_src (datatype *i_data_source, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : real_plan <datatype> (i_data_in, i_data_out, i_coeff), data_source (i_data_source) {
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public real_plan <datatype>::factory
			{
			private:
				datatype *data_source; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (datatype *i_data_source, datatype i_coeff = 1.0) : real_plan <datatype>::factory (i_coeff), data_source (i_data_source) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new z_src <datatype> (data_source, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a source term to an equation
		 ************************************************************************/
		template <class datatype>
		class pressure_grad_1d : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_out;
			using explicit_plan <datatype>::grid_m;
		
			const datatype *data_grad; //!< The data pointer for the source data
			datatype *data_top; //!< The data pointer for the source data
			datatype *data_bot; //!< The data pointer for the source data
			const datatype *pos_m;
			std::vector <datatype> diff;
		
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			pressure_grad_1d (grids::variable <datatype> &i_data_top, grids::variable <datatype> &i_data_bot, datatype *i_data_grad, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
			explicit_plan <datatype> (i_data_in, i_data_out, i_coeff), 
			data_grad (i_data_grad),
			data_top (i_data_top.ptr ()),
			data_bot (i_data_bot.ptr ()) {
				TRACE ("Adding source...");
				pos_m = &grid_m [0];

				diff.resize (m);
				for (int j = 1; j < m - 1; ++j)
				{
					diff [j] = pos_m [j + 1] - pos_m [j - 1];
				}
				diff [0] = pos_m [1] - pos_m [0];
				diff [m - 1] = pos_m [m - 1] - pos_m [m - 2];
			}
		
			virtual ~pressure_grad_1d () {}

			/*!**********************************************************************
			 * \copydoc explicit_plan::execute
			 ************************************************************************/
			virtual void execute () {
				TRACE ("Executing source...");
				for (int i = 0; i < ldn; ++i)
				{
					for (int j = 1; j < m - 1; ++j)
					{
						data_out [i * m + j] += coeff * data_top [i * m + j] / data_bot [i * m + j] * (data_grad [j + 1] - data_grad [j - 1]) / diff [j];
					}
					data_out [i * m] += coeff * data_top [i * m] / data_bot [i * m] * (data_grad [1] - data_grad [0]) / diff [0];
					data_out [i * m + m - 1] += coeff * data_top [i * m + m - 1] / data_bot [i * m + m - 1] * (data_grad [m - 1] - data_grad [m - 2]) / diff [m - 1];
				}	
			}
		
			/*!**********************************************************************
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			private:
				datatype *data_grad; //!< The data source to be used when constructing the plan
				grids::variable <datatype> &data_top; //!< The data source to be used when constructing the plan
				grids::variable <datatype> &data_bot; //!< The data source to be used when constructing the plan
			
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (grids::variable <datatype> &i_data_top, grids::variable <datatype> &i_data_bot, datatype *i_data_grad, datatype i_coeff = 1.0) : 
				explicit_plan <datatype>::factory (i_coeff), 
				data_grad (i_data_grad),
				data_top (i_data_top),
				data_bot (i_data_bot) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new pressure_grad_1d <datatype> (data_top, data_bot, data_grad, i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};

		/*!**********************************************************************
		 * \brief A plan to add a constant term to an equation
		 ************************************************************************/
		template <class datatype>
		class constant : public explicit_plan <datatype>
		{
		private:
			using explicit_plan <datatype>::coeff;
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::dims;
			using explicit_plan <datatype>::data_out;
				
		public:
			/*!**********************************************************************
			 * \copydoc explicit_plan::explicit_plan
			 * 
			 * \param i_coeff The coefficient for the source term
			 * \param i_data_source The data pointer for the source data
			 * 
			 * In this plan, data_source is not used in leiu of data_in. The reason for this is that data_in is almost always assumed to be the current variable rather than some other source term.
			 ************************************************************************/
			constant (grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : explicit_plan <datatype> (i_data_in, i_data_out, i_coeff) {
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
			 * \copydoc explicit_plan::factory
			 ************************************************************************/
			class factory : public explicit_plan <datatype>::factory
			{
			public:
				/*!**********************************************************************
				 * \param i_coeff The coefficient to be used when constructing the plan
				 * \param i_data_source The data source to be used when constructing the plan
				 ************************************************************************/
				factory (datatype i_coeff = 1.0) : explicit_plan <datatype>::factory (i_coeff) {}
			
				virtual ~factory () {}
			
				/*!**********************************************************************
				 * \copydoc explicit_plan::factory::_instance
				 ************************************************************************/
				virtual std::shared_ptr <plans::plan <datatype> > _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					if (coeff) {
						return std::shared_ptr <plans::plan <datatype> > (new constant <datatype> (i_data_in, i_data_out, coeff));
					}
					return std::shared_ptr <plans::plan <datatype> > ();
				}
			};
		};
	} /* source */
} /* plans */

#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

