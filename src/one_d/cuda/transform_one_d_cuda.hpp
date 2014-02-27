/*!***********************************************************************
 * \file transform_one_d_cuda.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_ONE_D_CUDA_HPP_G5118SR0
#define FFTW_ONE_D_CUDA_HPP_G5118SR0

#include "../../config.hpp"
#include "../../bases/plan.hpp"
#include "../../utils/cuda/utils_cuda.hpp"
#include "../../bases/transform.hpp"

namespace cuda
{
	namespace one_d
	{
		/*!*******************************************************************
		 * An implementation of the transform class using FFTW3.
		 * 
		 * \brief \copybrief bases::transform
		 *********************************************************************/
		template <class datatype>
		class transform : public bases::plan <datatype>
		{
		public:
			/*!*******************************************************************
			 * \copydoc bases::transform::transform ()
			 *********************************************************************/
			transform (bases::grid <datatype> &i_grid, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags);

			virtual ~transform ();

			/*!*******************************************************************
			 * \copydoc bases::transform::execute ()
			 *********************************************************************/
			void execute ();

		protected:
			int n;
			void *data_in;
			void *data_complex;
			void *data_out;

			datatype scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
			void *cu_plan;
		};
	
		template <class datatype>
		class master_transform : public bases::master_transform <datatype>
		{
		public:
			master_transform (bases::grid <datatype> &i_grid, datatype* i_data_in, datatype* i_data_out, int i_flags, int *element_flags, int *component_flags):
			bases::master_transform <datatype> (element_flags, component_flags),
			ldn (i_grid.ld),
			data_in (i_data_in),
			data_out (i_data_out ? i_data_out : i_data_in) {
				data.resize (ldn);
				if (i_flags & forward_vertical) {
					forward_transform = new one_d::transform <datatype> (i_grid, &data [0], &data [0], 0x00, element_flags, component_flags);
				}
				if (i_flags & inverse_vertical) {
					if (forward_transform) {
						inverse_transform = forward_transform;
					} else {
						inverse_transform = new one_d::transform <datatype> (i_grid, &data [0], &data [0], inverse, element_flags, component_flags);
					}
				}
			}
		
			virtual ~master_transform () {
				delete (forward_transform);
				delete (inverse_transform);
			}
	
			void transform (int flags) {
				if (flags & forward_vertical) {
					if (!(*component_flags & transformed_vertical) && forward_transform) {
						forward_transform->execute ();
					}
				}
				if (flags & inverse_vertical) {
					if ((*component_flags & transformed_vertical) && inverse_transform) {
						inverse_transform->execute ();
					}
				}
			}
	
			void write () {
				data.copy_to_device (ldn, data_in);
			}

			void read () {
				data.copy_to_host (ldn, data_out);
			}
		
		private:
			int &ldn;
			datatype *data_in, *data_out;
			cuda::vector <datatype> data;
	
			bases::plan <datatype> *forward_transform;
			bases::plan <datatype> *inverse_transform;
	
			using bases::master_transform <datatype>::component_flags;
		};
	} /* one_d */
} /* cuda */

#endif /* end of include guard: FFTW_ONE_D_CUDA_HPP_G5118SR0 */
