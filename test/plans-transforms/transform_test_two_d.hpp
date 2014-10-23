/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>
#include <fftw3.h>

#define TEST_TINY 1.e-4

class transform_two_d_test_suite : public CxxTest::TestSuite
{
public:
	void test_transform_chebyshev_fourier () {
		int n = 100, m = 200, ldn = 2 * (n / 2 + 1);
		std::vector <double> init (ldn * m, 0.0), init_copy (ldn * m, 0.0), fftw_in (ldn * m, 0.0), fftw_out (ldn * m, 0.0);
		double scalar;
		std::string file_name = "output";
	
		scalar = 1.0 / std::sqrt (n * 2 * (m - 1));

		fftw_iodim inverse_iodim, inverse_major_iodim, forward_iodim, forward_major_iodim;
		fftw_plan forward_horizontal_transform, inverse_horizontal_transform, vertical_transform;
		
		forward_iodim.n = n;
		forward_major_iodim.is = 1;
		forward_major_iodim.os = 1;
		
		forward_iodim.is = m;
		forward_iodim.os = 2 * m;
		forward_major_iodim.n = m;
		forward_horizontal_transform = fftw_plan_guru_split_dft_r2c (1, &forward_iodim, 1, &forward_major_iodim, &fftw_in [0], &fftw_out [0], &fftw_out [0] + m, FFTW_MEASURE);
		
		inverse_iodim.n = n;
		inverse_major_iodim.is = 1;
		inverse_major_iodim.os = 1;

		inverse_iodim.is = 2 * m;
		inverse_iodim.os = m;
		inverse_major_iodim.n = m;
		inverse_horizontal_transform = fftw_plan_guru_split_dft_c2r (1, &inverse_iodim, 1, &inverse_major_iodim, &fftw_out [0], &fftw_out [0] + m, &fftw_in [0], FFTW_MEASURE);
		
		fftw_r2r_kind kind = FFTW_REDFT00;
		vertical_transform = fftw_plan_many_r2r (1, &m, ldn, &fftw_out [0], NULL, 1, m, &fftw_out [0], NULL, 1, m, &kind, FFTW_MEASURE);
		
		srand (1);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				init [i * m + j] = rand () % 100;
				init_copy [i * m + j] = init [i * m + j];
				fftw_in [i * m + j] = init [i * m + j];
			}
		}
		
		fftw_execute (forward_horizontal_transform);
		fftw_execute (vertical_transform);
		fftw_execute (vertical_transform);
		fftw_execute (inverse_horizontal_transform);
		
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("Transform failure in FFTW (This should never happen)", init_copy [i * m + j], fftw_in [i * m + j] * scalar * scalar, TEST_TINY);
			}
		}
	}
};
