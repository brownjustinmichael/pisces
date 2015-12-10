/*!**********************************************************************
 * \file plans_variable_test.hpp
 * /Users/justinbrown/Dropbox/pisces/test/plans/plans_variable_test.hpp
 ************************************************************************/

#include <cxxtest/TestSuite.h>

#include "plans/grids/variable.hpp"

class plans_variable_test_suite : public CxxTest::TestSuite
{
public:
	void test_multiplication () {
		int n = 200, m = 300, flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);
		
		grids::variable a (grid_n, grid_m, flags), b (grid_n, grid_m, flags), c (grid_n, grid_m, flags), d (grid_n, grid_m, flags);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				a [i * m + j] = i;
				b [i * m + j] = j;
				d [i * m + j] = i * j;
			}
		}

		c == a * b;

		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Multiplication failure", c [i * m + j], d [i * m + j], 1.0e-8);
			}
		}
	}

	void test_division () {
		int n = 200, m = 300, flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);
		
		grids::variable a (grid_n, grid_m, flags), b (grid_n, grid_m, flags), c (grid_n, grid_m, flags), d (grid_n, grid_m, flags);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				a [i * m + j] = i;
				b [i * m + j] = j;
				d [i * m + j] = ((double) i) / ((double) j);
			}
		}

		c == a / b;

		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Multiplication failure", c [i * m + j], d [i * m + j], 1.0e-8);
			}
		}
	}


	void test_addition () {
		int n = 200, m = 300, flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);
		
		grids::variable a (grid_n, grid_m, flags), b (grid_n, grid_m, flags), c (grid_n, grid_m, flags), d (grid_n, grid_m, flags);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				a [i * m + j] = i;
				b [i * m + j] = j;
				d [i * m + j] = i + j;
			}
		}

		c == a + b;

		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Multiplication failure", c [i * m + j], d [i * m + j], 1.0e-8);
			}
		}
	}

	void test_subtraction () {

		int n = 200, m = 300, flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);
		
		grids::variable a (grid_n, grid_m, flags), b (grid_n, grid_m, flags), c (grid_n, grid_m, flags), d (grid_n, grid_m, flags);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				a [i * m + j] = i;
				b [i * m + j] = j;
				d [i * m + j] = i - j;
			}
		}

		c == a - b;

		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Multiplication failure", c [i * m + j], d [i * m + j], 1.0e-8);
			}
		}
	}

	void test_pemdas () {
		int n = 200, m = 300, flags;

		grids::horizontal::grid grid_n (n, -1.0, 1.0);
		grids::vertical::grid grid_m (m, -1.0, 1.0);
		
		grids::variable a (grid_n, grid_m, flags), b (grid_n, grid_m, flags), c (grid_n, grid_m, flags), d (grid_n, grid_m, flags);

		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				a [i * m + j] = i;
				b [i * m + j] = j;
				d [i * m + j] = 1. / (i * j + j);
			}
		}

		c == 1. / (a * b + b);

		for (int i = 1; i < n - 1; ++i) {
			for (int j = 1; j < m - 1; ++j) {
				TSM_ASSERT_DELTA ("Multiplication failure", c [i * m + j], d [i * m + j], 1.0e-8);
			}
		}
	}
};
