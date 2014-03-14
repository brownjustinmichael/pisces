/*!**********************************************************************
 * \file transform_two_d_cuda.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_CUDA_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_CUDA_HPP_RAXYBFTC

#include "../config.hpp"
#include "plan_two_d.hpp"
#include "../bases/transform.hpp"
#include "../../utils/cuda/utils_cuda.hpp"
#include "../utils/utils.hpp"

namespace two_d
{
	namespace fourier
	{
		template <class datatype>
		class horizontal_transform : public bases::plan <datatype>
		{
		public:
			/*!**********************************************************************
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			horizontal_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags = 0);

			horizontal_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags = 0);
			
			virtual ~horizontal_transform () {}
			
			virtual void execute ();
		
		protected:
			int n;
			int m;
			void *data_in;
			void *data_out;
			
			int flags;
			int threads;
			datatype scalar;
			void *cu_plan;
		};
		
		template <class datatype>
		class vertical_transform : public bases::plan <datatype>
		{
		public:
			/*!**********************************************************************
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			vertical_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags = 0);

			vertical_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags = 0);
			
			virtual ~vertical_transform () {}
			
			virtual void execute ();
		
		protected:
			int n;
			int m;
			void *data_in;
			void *data_complex;
			void *data_out;
			
			int flags;
			int threads;
			datatype scalar;
			void *cu_plan;
		};
		
		template <class datatype>
		class master_transform : public bases::master_transform <datatype>
		{
		public:
			master_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags) : 
			bases::master_transform <datatype> (element_flags, component_flags),
			ldn (i_grid_n.ld),
			ldm (i_grid_m.ld),
			data_in (i_data_in),
			data_out (i_data_out ? i_data_out : i_data_in) {
				data.resize (ldn * ldm);
				if (i_flags & forward_vertical) {
					forward_vertical_transform = new two_d::fourier::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, component_flags, i_threads);
				}
				if (i_flags & inverse_vertical) {
					if (forward_vertical_transform) {
						inverse_vertical_transform = forward_vertical_transform;
					} else {
						inverse_vertical_transform = new two_d::fourier::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, component_flags, i_threads);
					}
				}
				if (i_flags & forward_horizontal) {
					forward_horizontal_transform = new two_d::fourier::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, component_flags, i_threads);
				}
				if (i_flags & inverse_horizontal) {
					inverse_horizontal_transform = new two_d::fourier::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, component_flags, i_threads);
				}
			}
			
			virtual ~master_transform () {
				delete (forward_horizontal_transform);
				delete (forward_vertical_transform);
				delete (inverse_horizontal_transform);
				delete (inverse_vertical_transform);
			}
		
			void _transform (int flags) {
				if (flags & forward_horizontal) {
					if (!(*component_flags & transformed_horizontal) && forward_horizontal_transform) {
						forward_horizontal_transform->execute ();
					}
				}
				if (flags & forward_vertical) {
					if (!(*component_flags & transformed_vertical) && forward_vertical_transform) {
						forward_vertical_transform->execute ();
					}
				}
				if (flags & inverse_horizontal) {
					if ((*component_flags & transformed_horizontal) && inverse_horizontal_transform) {
						inverse_horizontal_transform->execute ();
					}
				}
				if (flags & inverse_vertical) {
					if ((*component_flags & transformed_vertical) && inverse_vertical_transform) {
						inverse_vertical_transform->execute ();
					}
				}
			}
		
			void write () {
				utils::matrix_copy (ldm, ldn, data_in, &data [0]);
			}
	
			void read () {
				utils::matrix_copy (ldm, ldn, &data [0], data_out);
			}
		
		private:
			int ldn, ldm;
			datatype *data_in, *data_out;
			utils::cuda::vector <datatype> data;
			
			bases::plan <datatype> *forward_horizontal_transform;
			bases::plan <datatype> *forward_vertical_transform;
			bases::plan <datatype> *inverse_horizontal_transform;
			bases::plan <datatype> *inverse_vertical_transform;
			
			using bases::master_transform <datatype>::component_flags;
		};
	} /* fourier */
} /* two_d */

#endif /* end of include guard: TRANSFORM_TWO_D_CUDA_HPP_RAXYBFTC */
