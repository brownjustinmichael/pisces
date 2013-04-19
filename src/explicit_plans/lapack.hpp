/*!***********************************************************************
 * \file lapack.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LAPACK_HPP_NSEHAOMK
#define LAPACK_HPP_NSEHAOMK

/*!*******************************************************************
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_(int *n, double *x, int *incx, double *y, int *incy);

namespace explicit_plans
{
	class copy : public plan
	{
	public:
		copy (int i_n, double *i_data_in, double *i_data_out) {
			n = i_n;
			data_in = i_data_in;
			data_out = i_data_out;
		}
		virtual ~copy () {}
		inline virtual void execute () {
			int ione = 1;
			dcopy_ (&n, data_in, &ione, data_out, &ione);
		}
	
		inline static std::unique_ptr<plan> make_unique (int i_n, double *i_data_in, double *i_data_out) {
			return std::unique_ptr<plan> (new copy (i_n, i_data_in, i_data_out));
		}

	private:
		int n;
		double *data_in;
		double *data_out;
	};

	class zero : public plan
	{
	public:
		zero (int i_n, double *i_data) {
			n = i_n;
			data = i_data;
		}
		virtual ~zero () {}
		inline virtual void execute () {
			for (int i = 0; i < n; ++i) {
				data [i] = 0.0;
			}
		}
	
		inline static std::unique_ptr<plan> make_unique (int i_n, double * i_data) {
			return std::unique_ptr<plan> (new zero (i_n, i_data));
		}
	
	private:
		int n;
		double *data;
	};
} /* explicit_plans */

#endif /* end of include guard: LAPACK_HPP_NSEHAOMK */
