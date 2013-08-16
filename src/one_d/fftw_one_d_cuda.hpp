/*!***********************************************************************
 * \file fftw_one_d_cuda.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FFTW_ONE_D_CUDA_HPP_G5118SR0
#define FFTW_ONE_D_CUDA_HPP_G5118SR0

#include "../config.hpp"
#include "../bases/transform.hpp"

namespace bases
{
	class element;
} /* bases */

namespace one_d
{
	namespace cuda
	{
		/*!*******************************************************************
		 * An implementation of the transform class using FFTW3.
		 * 
		 * \brief \copybrief bases::transform
		 *********************************************************************/
		class fftw_cosine : public bases::transform
		{
		public:
			/*!*******************************************************************
			 * \copydoc bases::transform::transform ()
			 *********************************************************************/
			fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out = null);
		
			virtual ~fftw_cosine () {}
		
			/*!*******************************************************************
			 * \copydoc bases::transform::execute ()
			 *********************************************************************/
			void execute ();
	
		private:
			int padded_n;
			double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
			fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		};
	} /* cuda */
} /* one_d */

#endif /* end of include guard: FFTW_ONE_D_CUDA_HPP_G5118SR0 */
