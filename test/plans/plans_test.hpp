/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>
#include "plans/advection.hpp"

#define TEST_TINY 1.e-4

class plans_test_suite : public CxxTest::TestSuite
{
public:
	void test_advection () {
		int n = 100, m = 200, ldn = 2 * (n / 2 + 1);
		std::vector <double> init (ldn * m, 0.0), init_copy (ldn * m, 0.0);
		
		
		
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("Transform failure in FFTW (This should never happen)", init_copy [i * m + j], init [i * m + j], TEST_TINY);
			}
		}
	}
};
