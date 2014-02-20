/*!**********************************************************************
 * \file element_one_d_cuda.hpp
 * /Users/justinbrown/Dropbox/pisces/src
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
	} /* chebyshev */
} /* one_d */

namespace cuda
{
	namespace one_d
	{
		namespace chebyshev
		{
			template <class datatype>
			class fft_element : public ::one_d::chebyshev::element <datatype>
			{
			public:
				fft_element (int i_n, int i_excess_0, datatype i_position_0, int i_excess_n, datatype i_position_n, int i_name, io::parameters& i_input_Params, bases::messenger* i_messenger_ptr, int i_flags);
		
				virtual ~fft_element () {}
		
				virtual void setup ();
	
				inline void explicit_reset () {
					::one_d::chebyshev::element <datatype>::explicit_reset ();
					
					TRACE ("Resetting explicit part...");
					
					utils::scale (n, (datatype) 0.0, rhs_dev.ptr ());
				}
				
				inline void implicit_reset () {
					::one_d::chebyshev::element <datatype>::implicit_reset ();
			
					if (!(flags & factorized)) {
						::utils::scale (n * n, 0.0, &matrix [0]);
					}
				}
		
				virtual datatype calculate_timestep ();
	
			private:
				using ::one_d::chebyshev::element <datatype>::n;
				using ::one_d::chebyshev::element <datatype>::flags;
				using ::one_d::chebyshev::element <datatype>::name;
				using ::one_d::chebyshev::element <datatype>::normal_stream;
				using ::one_d::chebyshev::element <datatype>::transform_stream;
				using ::one_d::chebyshev::element <datatype>::cell;
				using ::one_d::chebyshev::element <datatype>::timestep;
				using ::one_d::chebyshev::element <datatype>::boundary_weights;
				using ::one_d::chebyshev::element <datatype>::params;
				using ::one_d::chebyshev::element <datatype>::grid;
				using ::one_d::chebyshev::element <datatype>::ptr;
				using ::one_d::chebyshev::element <datatype>::messenger_ptr;

				int excess_0, excess_n;
				utils::vector <datatype> data_dev;
				utils::vector <datatype> rhs_dev;
				std::vector<datatype> matrix; //!< A vector containing the datatype matrix used in the implicit solver
			};
		} /* chebyshev */
	} /* one_d */
} /* cuda */

#endif /* end of include guard: ELEMENT_ONE_D_CUDA_HPP_IBOLYZP9 */
