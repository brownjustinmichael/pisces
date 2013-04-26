/*!***********************************************************************
 * \file one_d/lapack.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LAPACK_HPP_NSEHAOMK
#define LAPACK_HPP_NSEHAOMK

#include "../bases/plan.hpp"

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_ (int *n, double *dx, int *incx, double *dy, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that scales a double array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param da The double by which to scale the data
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void dscal_ (int *n, double *da, double *dx, int *incx);

namespace one_d
{
	/*!*******************************************************************
	 * \brief A plan that scales an array into another or in place
	 *********************************************************************/
	class scale : public bases::explicit_plan
	{
	public:	
		/*!*******************************************************************
		 * \param i_scalar The double by which to scale the array
		 * \copydoc bases::explicit_plan::explicit_plan ()
		 *********************************************************************/
		scale (double i_scalar, int i_n, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr, i_logger) {
			scalar = i_scalar;
		}
		
		scale (double i_scalar, int i_n, double &i_data_in, double &i_data_out, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, &i_data_out, i_flags_ptr, i_logger) {
			scalar = i_scalar;
		}
		
		scale (double i_scalar, int i_n, double &i_data_in, int *i_flags_ptr = NULL, int i_logger = -1) : bases::explicit_plan (i_n, &i_data_in, NULL, i_flags_ptr, i_logger) {
			scalar = i_scalar;
		}
		
		virtual ~scale () {}
		
		/*!*******************************************************************
		 * \copydoc bases::explicit_plan::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			int ione = 1;
			if (data_out != data_in) {
				dcopy_ (&n, data_in, &ione, data_out, &ione);
			}
			
			dscal_ (&n, &scalar, data_in, &ione);
		}

	private:
		double scalar; //!< The double by which to scale the data
	};
} /* one_d */

#endif /* end of include guard: LAPACK_HPP_NSEHAOMK */
