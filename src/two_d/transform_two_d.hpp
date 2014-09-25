/*!**********************************************************************
 * \file transform_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-29.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef TRANSFORM_TWO_D_HPP_RAXYBFTC
#define TRANSFORM_TWO_D_HPP_RAXYBFTC

#include "logger/logger.hpp"
#include "plan_two_d.hpp"
#include "../bases/transform.hpp"
#include <fftw3.h>
#include "linalg/utils.hpp"

namespace two_d
{
	// class fftw_configure
	// {
	// public:
	// 	fftw_configure ();
	// 	
	// 	virtual ~fftw_configure ();
	// 	
	// private:
	// 	bool threads;
	// };
	
	namespace fourier
	{
		template <class datatype>
		class horizontal_transform : public bases::plan <datatype>
		{
		public:
			/*!**********************************************************************
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			horizontal_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);

			horizontal_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);
			
			virtual ~horizontal_transform () {}
			
			virtual void execute ();
		
		protected:
			int n;
			int m;
			datatype *data_in;
			datatype *data_out;
			
			int flags;
			int threads;
			datatype scalar;
			std::vector <fftw_plan> plans;
			std::vector <fftwf_plan> plans_float;
			fftw_iodim major_iodim;
			fftw_iodim iodim;
		};
		
		template <class datatype>
		class vertical_transform : public bases::plan <datatype>
		{
		public:
			/*!**********************************************************************
			 * WARNING!! BECAUSE OF THE REAL DATA FFT, THE ARRAYS MUST HAVE DIMENSION M * 2 * (N / 2 + 1)
			 ************************************************************************/
			vertical_transform (int n, int m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);

			vertical_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *i_element_flags, int *i_component_flags, int i_threads = 0);
			
			virtual ~vertical_transform () {}
			
			virtual void execute ();
		
		protected:
			int n;
			int m;
			datatype *data_in;
			datatype *data_out;
			
			int flags;
			int threads;
			datatype scalar;
			std::vector <fftw_plan> plans;
			std::vector <fftwf_plan> plans_float;
		};
		
		template <class datatype>
		class master_transform : public bases::master_transform <datatype>
		{
		public:
			master_transform (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags, int i_threads) : 
			bases::master_transform <datatype> (element_flags, component_flags),
			ldn (i_grid_n.get_ld ()),
			ldm (i_grid_m.get_ld ()),
			data_in (i_data_in),
			data_out (i_data_out ? i_data_out : i_data_in) {
				data.resize (ldn * ldm, 0.0);
				if (i_flags & forward_vertical) {
					forward_vertical_transform = std::shared_ptr <bases::plan <datatype>> (new two_d::fourier::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, component_flags, i_threads));
				}
				if (i_flags & inverse_vertical) {
					if (forward_vertical_transform) {
						inverse_vertical_transform = forward_vertical_transform;
					} else {
						inverse_vertical_transform = std::shared_ptr <bases::plan <datatype>> (new two_d::fourier::vertical_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, component_flags, i_threads));
					}
				}
				if (i_flags & forward_horizontal) {
					forward_horizontal_transform = std::shared_ptr <bases::plan <datatype>> (new two_d::fourier::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], 0x00, element_flags, component_flags, i_threads));
				}
				if (i_flags & inverse_horizontal) {
					inverse_horizontal_transform = std::shared_ptr <bases::plan <datatype>> (new two_d::fourier::horizontal_transform <datatype> (i_grid_n, i_grid_m, &data [0], &data [0], inverse, element_flags, component_flags, i_threads));
				}
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
				// std::stringstream debug;
				// for (int i = 0; i < 2 * (32 / 2 + 1); ++i) {
				// 	debug << data_in [i * 32 + 32 - 1] << " ";
				// }
				// DEBUG ("WRITING " << debug.str ());
				// debug.str ("");
				utils::matrix_copy (ldm, ldn, data_in, &data [0]);
			}
	
			void read () {
				utils::matrix_copy (ldm, ldn, &data [0], data_out);
				// std::stringstream debug;
				// for (int i = 0; i < 2 * (32 / 2 + 1); ++i) {
				// 	debug << data_out [i * 32 + 32 - 1] << " ";
				// }
				// DEBUG ("READING " << debug.str ());
				// debug.str ("");
			}
			
			virtual datatype *get_data () {
				return &data [0];
			}
		
		private:
			int ldn, ldm;
			datatype *data_in, *data_out;
			std::vector <datatype> data;
			
			std::shared_ptr <bases::plan <datatype>> forward_horizontal_transform;
			std::shared_ptr <bases::plan <datatype>> forward_vertical_transform;
			std::shared_ptr <bases::plan <datatype>> inverse_horizontal_transform;
			std::shared_ptr <bases::plan <datatype>> inverse_vertical_transform;
			
			using bases::master_transform <datatype>::component_flags;
		};
	} /* fourier */
} /* two_d */

#endif /* end of include guard: TRANSFORM_TWO_D_HPP_RAXYBFTC */
