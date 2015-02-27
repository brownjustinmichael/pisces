/*!**********************************************************************
 * \file block_test.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-11-05.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include "mpi/messenger.hpp"
#include <vector>
#include <stdio.h>
#include "linalg/linalg.hpp"
#include "plans/grid.hpp"
#include "plans/advection.hpp"

#define TEST_TINY 1.e-6

class plan_solver_test_suite : public CxxTest::TestSuite
{
private:
	mpi::messenger mess;
	
public:
	void test_advection () {
		std::vector <double> x_velocity, z_velocity, data, rhs, rhs_compare;
		int n = 200, m = 500;
		
		plans::vertical::grid ()
	}
};

