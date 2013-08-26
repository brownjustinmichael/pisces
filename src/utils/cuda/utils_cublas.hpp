/*!**********************************************************************
 * \file utils_cublas.hcu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_CUBLAS_HCU_K6B7R3SS
#define UTILS_CUBLAS_HCU_K6B7R3SS

namespace cuda
{
	namespace utils
	{
		template <class datatype>
		struct vector 
		{
		public:
			vector (int i_n = 0, datatype *x = NULL, int incx = 1);
		
			~vector ();
		
			datatype* pointer () {
				return vect;
			}
			
			void resize (int i_n);
		
			int size () {
				return n;
			}
		
			void copy_to_device (int n, datatype* x, int incx = 1);
		
			void copy_to_host (int n, datatype* x, int incx = 1);
		
		private:
			int n;
			datatype* vect;
			int inc_vect;
		};
		
		void copy (int n, float* x, float* y, int incx = 1, int incy = 1);
	
		void copy (int n, double* x, double* y, int incx = 1, int incy = 1);
	
		void scale (int n, float a, float* x, int incx = 1);
		
		void scale (int n, double a, double* x, int incx = 1);

		double dot (int n, float* x, float* y, int incx = 1, int incy = 1);
	
		double dot (int n, double* x, double* y, int incx = 1, int incy = 1);

		void add_scaled (int n, float a, float* x, float* y, int incx = 1, int incy = 1);
	
		void add_scaled (int n, double a, double* x, double* y, int incx = 1, int incy = 1);
	
		void matrix_vector_multiply (int m, int n, float alpha, float* a, float* x, float beta, float* y, int lda, int incx = 1, int incy = 1);
		
		void matrix_vector_multiply (int m, int n, double alpha, double* a, double* x, double beta, double y, int lda, int incx = 1, int incy = 1);
	} /* utils */
} /* cuda */

#endif /* end of include guard: UTILS_CUBLAS_HCU_K6B7R3SS */
