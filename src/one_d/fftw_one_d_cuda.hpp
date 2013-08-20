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
#include "../bases/plan.hpp"

#ifndef __CUDACC__
class cufftHandle;
class cufftDoubleReal;
class cufftDoubleComplex;
#endif // __CUDACC__

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
		 * \brief \copybrief bases::explicit_plan
		 *********************************************************************/
		class fftw_cosine : public bases::explicit_plan
		{
		public:
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan::explicit_plan ()
			 *********************************************************************/
			fftw_cosine (bases::element* i_element_ptr, int i_n, int i_name_in, int i_name_out = null);
		
			virtual ~fftw_cosine ();
		
			/*!*******************************************************************
			 * \copydoc bases::explicit_plan::execute ()
			 *********************************************************************/
			void execute ();
	
		private:
			double scalar; //!< The scalar used after the transform (1 / sqrt (2 * (n - 1)))
			cufftDoubleReal* data_real;
			cufftDoubleComplex* data_complex;
			cufftHandle* cu_plan;
			// fftw_plan fourier_plan; //!< The fftw_plan object to be executed
		};
	} /* cuda */
} /* one_d */

#endif /* end of include guard: FFTW_ONE_D_CUDA_HPP_G5118SR0 */
