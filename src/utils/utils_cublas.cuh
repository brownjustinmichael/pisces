/*!**********************************************************************
 * \file utils_cublas.hcu
 * /Users/justinbrown/Dropbox/spectral_element
 * 
 * Created by Justin Brown on 2013-08-15.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef UTILS_CUBLAS_HCU_K6B7R3SS
#define UTILS_CUBLAS_HCU_K6B7R3SS

#include <cublas.h>
#include <stdio.h>

namespace utils
{
	namespace cuda
	{
		void HANDLE_STATUS (cublasStatus status);
		
		struct config
		{
		public:
			config () {
				HANDLE_STATUS (cublasInit ());
			}
		
			virtual ~config () {
				HANDLE_STATUS (cublasShutdown ());
			}
		};
		
		template <class datatype>
		struct vector 
		{
		public:
			vector (int i_n, datatype *x = NULL, int incx = 1) {
				inc_vect = 1;
				n = i_n;
				HANDLE_STATUS (cublasAlloc (n, sizeof (datatype), (void**) &vect));
				if (x) {
					copy_to_device (n, x, incx);
				}
			}
		
			~vector () {
				HANDLE_STATUS (cublasFree (vect));
			}
		
			datatype* pointer () {
				return vect;
			}
		
			int size () {
				return n;
			}
		
			void copy_to_device (int n, datatype* x, int incx = 1) {
				HANDLE_STATUS (cublasSetVector (n, sizeof (datatype), x, incx, vect, inc_vect));
			}
		
			void copy_to_host (int n, datatype* x, int incx = 1) {
				HANDLE_STATUS (cublasGetVector (n, sizeof (datatype), vect, inc_vect, x, incx));
			}
		
		private:
			int n;
			datatype* vect;
			int inc_vect;
		};
		
		typedef vector <float> vect_float;
		typedef vector <double> vect_double;
		
		void copy (int n, vect_float* x, vect_float* y, int incx = 1, int incy = 1);
	
		void copy (int n, vect_double* x, vect_double* y, int incx = 1, int incy = 1);
	
		void scale (int n, float a, vect_float* x, int incx = 1);
		
		void scale (int n, double a, vect_double* x, int incx = 1);

		double dot (int n, vect_float* x, vect_float* y, int incx = 1, int incy = 1);
	
		double dot (int n, vect_double* x, vect_double* y, int incx = 1, int incy = 1);

		void add_scaled (int n, float a, vect_float* x, vect_float* y, int incx = 1, int incy = 1);
	
		void add_scaled (int n, double a, vect_double* x, vect_double* y, int incx = 1, int incy = 1);
	
		void matrix_vector_multiply (int m, int n, float alpha, vect_float* a, vect_float* x, float beta, vect_float* y, int lda, int incx = 1, int incy = 1);
		
		void matrix_vector_multiply (int m, int n, double alpha, vect_double* a, vect_double* x, double beta, vect_double y, int lda, int incx = 1, int incy = 1);
		
	} /* cublas */

} /* utils */

#endif /* end of include guard: UTILS_CUBLAS_HCU_K6B7R3SS */
