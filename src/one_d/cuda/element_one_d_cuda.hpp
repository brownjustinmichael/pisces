/*!**********************************************************************
 * \file element_one_d_cuda.hpp
 * /Users/justinbrown/Dropbox/spectral_element/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_ONE_D_CUDA_HPP_IBOLYZP9
#define ELEMENT_ONE_D_CUDA_HPP_IBOLYZP9

#include "../element_one_d.hpp"
#include "../../utils/cuda/utils_cublas.hpp"

namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		class element;
		
		namespace cuda
		{
			template <class datatype>
			class fft_element : public element <datatype>
			{
			public:
				fft_element (int i_n, datatype i_position_0, datatype i_position_n, int i_excess_0, int i_excess_n, int i_name, io::parameter_map& i_input_Params, bases::messenger <datatype>* i_messenger_ptr, int i_flags);
			
				virtual ~fft_element () {}
			
				virtual void setup ();
		
				inline void implicit_reset () {
					element <datatype>::implicit_reset ();
				
					if (!(flags & factorized)) {
						utils::scale (n * n, 0.0, &matrix [0]);
					}
				}
			
				virtual datatype calculate_timestep ();
		
			private:
				using element <datatype>::n;
				using element <datatype>::flags;
				using element <datatype>::name;
				using element <datatype>::normal_stream;
				using element <datatype>::transform_stream;
				using element <datatype>::cell;
				using element <datatype>::timestep;
				using element <datatype>::boundary_weights;
				using element <datatype>::inputParams;
				using element <datatype>::grid;

				int excess_0, excess_n;
				utils::cuda::vector <datatype> data_dev;
				utils::cuda::vector <datatype> rhs_dev;
				std::vector<datatype> matrix; //!< A vector containing the datatype matrix used in the implicit solver
			};
		} /* cuda */
	} /* chebyshev */
} /* one_d */

#endif /* end of include guard: ELEMENT_ONE_D_CUDA_HPP_IBOLYZP9 */
