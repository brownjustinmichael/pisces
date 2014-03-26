/*!**********************************************************************
 * \file source_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src/two_d
 * 
 * Created by Justin Brown on 2013-12-24.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef SOURCE_TWO_D_HPP_G9AN5CH6
#define SOURCE_TWO_D_HPP_G9AN5CH6

#include "plan_two_d.hpp"
#include <omp.h>

namespace two_d
{
	namespace fourier
	{
		template <class datatype>
		class source : public explicit_plan <datatype>
		{
		public:
			source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source) :
			explicit_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source (i_data_source) {}
			
			virtual ~source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");
				for (int j = 0; j < m; ++j) {
					utils::add_scaled (ldn, 1.0, data_source + j, data_out + j, m, m);
				}
			}
		
		private:
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source;
		};

		template <class datatype>
		class x_derivative_source : public explicit_plan <datatype>
		{
		public:
			x_derivative_source (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out, datatype i_coeff, datatype* i_data_source) :
			explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			coeff (i_coeff),
			data_source (i_data_source) {}
			
			x_derivative_source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source) :
			explicit_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source (i_data_source) {}
			
			virtual ~x_derivative_source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");
				std::stringstream debug;
				for (int i = 2; i < ldn; i += 2) {
					utils::add_scaled (m, coeff * 2.0 * acos (-1.0) / (grid_n [n - 1] - grid_n [0]) * (i / 2), data_source + i * m, data_out + (i + 1) * m);
					utils::add_scaled (m, coeff * -2.0 * acos (-1.0) / (grid_n [n - 1] - grid_n [0]) * (i / 2), data_source + (i + 1) * m, data_out + i * m);
				}
			}
		
		private:
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::grid_n;
			using explicit_plan <datatype>::grid_m;
			using explicit_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source;
		};
		
		template <class datatype>
		class square_x_derivative_source : public real_plan <datatype>
		{
		public:
			square_x_derivative_source (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out, datatype i_coeff, datatype* i_data_source) :
			real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_n (&(grid_n [0])) {}
			
			square_x_derivative_source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source) :
			real_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_n (&(grid_n [0])) {}
			
			virtual ~square_x_derivative_source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");
#pragma omp parallel for
				for (int j = 0; j < m; ++j) {
					data_out [j] += coeff * (data_source [m + j] - data_source [(n - 1) * m + j]) * (data_source [m + j] - data_source [(n - 1) * m + j]) / (pos_n [1] - pos_n [0]) / (pos_n [1] - pos_n [0]) / 4.0;
					for (int i = 1; i < n - 1; ++i) {
						data_out [i * m + j] += coeff * (data_source [(i + 1) * m + j] - data_source [(i - 1) * m + j]) * (data_source [(i + 1) * m + j] - data_source [(i - 1) * m + j]) / (pos_n [i + 1] - pos_n [i - 1]) / (pos_n [i + 1] - pos_n [i - 1]);
					}
					data_out [(n - 1) * m + j] += coeff * (data_source [j] - data_source [(n - 2) * m + j]) * (data_source [j] - data_source [(n - 2) * m + j]) / (pos_n [n - 1] - pos_n [n - 2]) / (pos_n [n - 1] - pos_n [n - 2]) / 4.0;
				}
				TRACE ("Execution complete.");
			}
		
		private:
			using real_plan <datatype>::n;
			using real_plan <datatype>::m;
			using real_plan <datatype>::grid_n;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source, *pos_n;
		};
		
		template <class datatype>
		class z_derivative_source : public explicit_plan <datatype>
		{
		public:
			z_derivative_source (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out, datatype i_coeff, datatype* i_data_source) :
			explicit_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_m (&(grid_m [0])) {}

			z_derivative_source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source) :
			explicit_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_m (&(grid_m [0])) {}
			
			virtual ~z_derivative_source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");
#pragma omp parallel for
				for (int i = 0; i < ldn; ++i) {
					data_out [i * m] += coeff * (data_source [i * m + 1] - data_source [i * m]) / (pos_m [1] - pos_m [0]);
					for (int j = 1; j < m - 1; ++j) {
						data_out [i * m + j] += coeff * (data_source [i * m + j + 1] - data_source [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
					}
					data_out [(i + 1) * m - 1] += coeff * (data_source [(i + 1) * m - 1] - data_source [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]);
				}
			}
		
		private:
			using explicit_plan <datatype>::n;
			using explicit_plan <datatype>::ldn;
			using explicit_plan <datatype>::m;
			using explicit_plan <datatype>::grid_m;
			using explicit_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source, *pos_m;
		};
		
		template <class datatype>
		class square_z_derivative_source : public real_plan <datatype>
		{
		public:
			square_z_derivative_source (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype *i_data_in, datatype *i_data_out, datatype i_coeff, datatype* i_data_source) :
			real_plan <datatype> (i_grid_n, i_grid_m, i_data_in, i_data_out),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_m (&(grid_m [0])) {}

			square_z_derivative_source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source) :
			real_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source (i_data_source),
			pos_m (&(grid_m [0])) {}
			
			virtual ~square_z_derivative_source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");
#pragma omp parallel for
				for (int i = 0; i < n; ++i) {
					data_out [i * m] += coeff * (data_source [i * m + 1] - data_source [i * m]) / (pos_m [1] - pos_m [0]) * (data_source [i * m + 1] - data_source [i * m]) / (pos_m [1] - pos_m [0]);
					for (int j = 1; j < m - 1; ++j) {
						data_out [i * m + j] += coeff * (data_source [i * m + j + 1] - data_source [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) * (data_source [i * m + j + 1] - data_source [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
					}
					data_out [(i + 1) * m - 1] += coeff * (data_source [(i + 1) * m - 1] - data_source [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) * (data_source [(i + 1) * m - 1] - data_source [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]);
				}
			}
		
		private:
			using real_plan <datatype>::n;
			using real_plan <datatype>::m;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source, *pos_m;
		};
		
		template <class datatype>
		class mixed_derivative_source : public real_plan <datatype>
		{
		public:
			mixed_derivative_source (bases::solver <datatype> &i_solver, datatype i_coeff, datatype* i_data_source_x, datatype* i_data_source_z) :
			real_plan <datatype> (i_solver),
			coeff (i_coeff),
			data_source_x (i_data_source_x),
			data_source_z (i_data_source_z),
			pos_n (&(grid_n [0])),
			pos_m (&(grid_m [0])) {}
			
			virtual ~mixed_derivative_source () {}
			
			virtual void execute () {
				TRACE ("Executing source...");

				data_out [0] += coeff * (data_source_z [1] - data_source_z [0]) / (pos_m [1] - pos_m [0]) * (data_source_x [m] - data_source_x [(n - 1) * m]) / (pos_n [1] - pos_n [0]) / 2.0;
				for (int j = 1; j < m - 1; ++j) {
					data_out [j] += coeff * (data_source_z [j + 1] - data_source_z [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) * (data_source_x [m + j] - data_source_x [(n - 1) * m + j]) / (pos_n [1] - pos_n [0]) / 2.0;
				}
				data_out [m - 1] += coeff * (data_source_z [m - 1] - data_source_z [m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) * (data_source_x [m - 1] - data_source_x [n * m - 1]) / (pos_n [1] - pos_n [0]) / 2.0;
#pragma omp parallel for
				for (int i = 1; i < n - 1; ++i) {
					data_out [i * m] += coeff * (data_source_z [i * m + 1] - data_source_z [i * m]) / (pos_m [1] - pos_m [0]) * (data_source_x [(i + 1) * m] - data_source_x [(i - 1) * m]) / (pos_n [i + 1] - pos_n [i - 1]);
					for (int j = 1; j < m - 1; ++j) {
						data_out [i * m + j] += coeff * (data_source_z [i * m + j + 1] - data_source_z [i * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) * (data_source_x [(i + 1) * m + j] - data_source_x [(i - 1) * m + j]) / (pos_n [i + 1] - pos_n [i - 1]);
					}
					data_out [(i + 1) * m - 1] += coeff * (data_source_z [(i + 1) * m - 1] - data_source_z [(i + 1) * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) * (data_source_x [(i + 2) * m - 1] - data_source_x [i * m - 1]) / (pos_n [i + 1] - pos_n [i - 1]);
				}
				data_out [(n - 1) * m] += coeff * (data_source_z [(n - 1) * m + 1] - data_source_z [(n - 1) * m]) / (pos_m [1] - pos_m [0]) * (data_source_x [0] - data_source_x [(n - 2) * m]) / (pos_n [n - 1] - pos_n [n - 2]) / 2.0;
				for (int j = 1; j < m - 1; ++j) {
					data_out [(n - 1) * m + j] += coeff * (data_source_z [(n - 1) * m + j + 1] - data_source_z [(n - 1) * m + j - 1]) / (pos_m [j + 1] - pos_m [j - 1]) * (data_source_x [j] - data_source_x [(n - 1) * m + j]) / (pos_n [n - 1] - pos_n [n - 2]) / 2.0;
				}
				data_out [n * m - 1] += coeff * (data_source_z [n * m - 1] - data_source_z [n * m - 2]) / (pos_m [m - 1] - pos_m [m - 2]) * (data_source_x [m - 1] - data_source_x [(n - 1) * m - 1]) / (pos_n [n - 1] - pos_n [n - 2]) / 2.0;
			}
		
		private:
			using real_plan <datatype>::n;
			using real_plan <datatype>::m;
			using real_plan <datatype>::grid_n;
			using real_plan <datatype>::grid_m;
			using real_plan <datatype>::data_out;
			datatype coeff;
			datatype *data_source_x, *data_source_z, *pos_n, *pos_m;
		};
	} /* fourier */
} /* two_d */


#endif /* end of include guard: SOURCE_TWO_D_HPP_G9AN5CH6 */

