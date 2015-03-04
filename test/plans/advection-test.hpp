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
		
		x_velocity.resize (n * m);
		z_velocity.resize (n * m);
		data.resize (n * m);
		rhs.resize (n * m);
		rhs_compare.resize (n * m);
		
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < n; ++i) {
				x_velocity [i * m + j] = rand () % 1000 / 1000;
				z_velocity [i * m + j] = rand () % 1000 / 1000;
				data [i * m + j] = rand () % 1000 / 1000;
			}
		}
		
		plans::horizontal::grid <double> grid_n (n, -1.0, 1.0);
		plans::vertical::grid <double> grid_m (m, -1.0, 1.0);
		
		plans::advection <double> plan (grid_n, grid_m, 1.0, &x_velocity [0], &z_velocity [0], &data [0], &rhs [0]);
		
		plan.execute ();
		
		for (int j = 1; j < m - 1; ++j) {
			for (int i = 1; i < n - 1; ++i) {
				rhs_compare [i * m + j] = x_velocity [i * m + j] * (data [(i + 1) * m + j] - data [(i - 1) * m + j]) / (grid_n [i + 1] - grid_n [i - 1]) + z_velocity [i * m + j] * (data [(i) * m + j + 1] - data [(i) * m + j - 1]) / (grid_m [j + 1] - grid_m [j - 1]);
			}
		}
		
		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Advection failure", rhs_compare [i * m + j], -rhs [i * m + j], 1.0e-8);
			}
		}
	}
};

