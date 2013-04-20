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
extern "C" void dcopy_(int *n, double *dx, int *incx, double *dy, int *incy);

namespace one_d
{
	/*!*******************************************************************
	 * \brief A plan that copies one array into another
	 *********************************************************************/
	class copy : public bases::explicit_plan
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::explicit_plan::explicit_plan ()
		 *********************************************************************/
		copy (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) : bases::explicit_plan (i_n, i_data_in, i_data_out, i_flags_ptr) {}
		
		virtual ~copy () {}
			
		/*!*******************************************************************
		 * \copydoc bases::explicit_plan::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			int ione = 1;
			dcopy_ (&n, data_in, &ione, data_out, &ione);
		}

		/*!*******************************************************************
		 * \brief Make a unique pointer to a new copy object
		 * \copydetails copy ()
		 *********************************************************************/
		inline static std::unique_ptr<plan> make_unique (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) {
			return std::unique_ptr<plan> (new copy (i_n, i_data_in, i_data_out, i_flags_ptr));
		}
	};

	/*!*******************************************************************
	 * \brief A plan that zeroes an array
	 *********************************************************************/
	class zero : public bases::plan
	{
	public:	
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data The double array of data
		 *********************************************************************/
		zero (int i_n, double *i_data) {
			n = i_n;
			data = i_data;
		}
		
		virtual ~zero () {}
		
		/*!*******************************************************************
		 * \copydoc bases::plan::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			for (int i = 0; i < n; ++i) {
				data [i] = 0.0;
			}
		}
		
		/*!*******************************************************************
		 * \brief Make a unique pointer to a new zero object
		 * \copydetails zero ()
		 *********************************************************************/
		inline static std::unique_ptr<plan> make_unique (int i_n, double * i_data) {
			return std::unique_ptr<plan> (new zero (i_n, i_data));
		}

	private:
		int n; //!< The integer number of elements in the data
		double *data; //!< The double array of data
	};
} /* one_d */

#endif /* end of include guard: LAPACK_HPP_NSEHAOMK */
