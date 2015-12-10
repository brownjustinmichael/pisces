/*!***********************************************************************
 * \file grid.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_CPP_HV4P0UOP
#define COLLOCATION_CPP_HV4P0UOP

#include <cassert>
#include <cmath>
#include "logger/logger.hpp"
#include "grid.hpp"

namespace grids
{
#ifndef _VCOS
	namespace vertical
	{
		grid::grid (axis *i_axis_ptr) : 
		grids::grid (i_axis_ptr, 3) {
			TRACE ("Instantiating...");

			scale = sqrt (2.0 / (n - 1));
			pioN = std::acos (-1.0) / (n - 1);
			exists_array.resize (n * n * 3, false);	
		
			if (n > 1) {
				double scale = (position_0 - position_n) / (std::cos (excess_0 * pioN) - std::cos ((n - 1 - excess_n) * pioN));
				double initial_position = position_0 - scale * std::cos (excess_0 * pioN);
				for (int i = 0; i < n; ++i) {
					positions [i] = scale * std::cos (i * pioN) + initial_position;
				}
			} else {
				positions [0] = (position_0 + position_n) / 2.0;
			}

			this->ood [0] = 1.0 / (positions [1] - positions [0]);
			this->ood2 [0] = 0.5 * this->ood [0];
			for (int i = 1; i < n - 1; ++i)
			{
				this->ood [i] = 1.0 / (positions [i + 1] - positions [i]);
				this->ood2 [i] = 1.0 / (positions [i + 1] - positions [i - 1]);
			}
			this->ood [n - 1] = 1.0 / (positions [n - 1] - positions [n - 2]);
			this->ood2 [n - 1] = 0.5 * this->ood [n - 1];
		
		
			width = positions [n - 1] - positions [0];

			// _calculate_matrix ();
	
			TRACE ("Instantiated...");
		}
		
		grid::grid (int i_n, double i_position_0, double i_position_n, int i_excess_0, int i_excess_n, int i_ld) : 
		grids::grid (3, i_n, i_position_0, i_position_n, i_excess_0, i_excess_n, i_ld) {
			TRACE ("Instantiating...");

			scale = sqrt (2.0 / (n - 1));
			pioN = std::acos (-1.0) / (n - 1);
			exists_array.resize (n * n * 3, false);	
		
			if (n > 1) {
				double scale = (position_0 - position_n) / (std::cos (excess_0 * pioN) - std::cos ((n - 1 - excess_n) * pioN));
				double initial_position = position_0 - scale * std::cos (excess_0 * pioN);
				for (int i = 0; i < n; ++i) {
					positions [i] = scale * std::cos (i * pioN) + initial_position;
				}
			} else {
				positions [0] = (position_0 + position_n) / 2.0;
			}
		
			this->ood [0] = 1.0 / (positions [1] - positions [0]);
			this->ood2 [0] = 0.5 * this->ood [0];
			for (int i = 1; i < n - 1; ++i)
			{
				this->ood [i] = 1.0 / (positions [i + 1] - positions [i]);
				this->ood2 [i] = 1.0 / (positions [i + 1] - positions [i - 1]);
			}
			this->ood [n - 1] = 1.0 / (positions [n - 1] - positions [n - 2]);
			this->ood2 [n - 1] = 0.5 * this->ood [n - 1];
			
		
			width = positions [n - 1] - positions [0];

			// _calculate_matrix ();
	
			TRACE ("Instantiated...");
		}

		double grid::recursion (int d, int m, int k) {		
			assert (d >= 0);
			assert (d < 3);
			assert (m >= 0);
			assert (k >= 0);

			// Use the recursion relations to calculate the correct value of the polynomial at the collocation point
			if (exists (d, m, k)) {
				// This value has already been calculated
				return grids::grid::index (d, m, k);
			} else if ((d == 2 && (m == 0 || m == 1)) || (d == 1 && m == 0)) {
				// The first two polynomials of the second derivative and the first polynomial of the first derivative are 0.0
				return 0.0;
			} else if ((d == 1 && m == 1) || (d == 0 && m == 0)) {
				// The second polynomial of the first derivative and the first Chebyshev polynomial are 1.0
				return scale;
			} else if (d == 0 && m == 1) {
				// The second Chebyshev polynomial is cos (theta)
				return scale * std::cos (pioN * k);
			} else if (d == 0) {
				// Use recursion to find the Chebyshev polynomial
				return 2.0 * std::cos (pioN * k) * recursion (0, m - 1, k) - recursion (0, m - 2, k);
			} else if (d == 1) {
				// Use recursion to find the first derivative of the Chebyshev polynomial
				return 2.0 * recursion (0, m - 1, k) + 2.0 * std::cos (pioN * k) * recursion (1, m - 1, k) - recursion (1, m - 2, k);
			} else if (d == 2) {
				// Use recursion to find the second derivative of the Chebyshev polynomial
				return 4.0 * recursion (1, m - 1, k) + 2.0 * std::cos (pioN * k) * recursion (2, m - 1, k) - recursion (2, m - 2, k);
			} else {
				// There's been a bad index, kill the program
				FATAL ("Bad index")
				throw 0;
			}
		}
		
		void grid::_calculate_matrix () {
			scale = sqrt (2.0 / (n - 1));
			pioN = std::acos (-1.0) / (n - 1);
			exists_array.resize (n * n * 3, false);	
		
			width = positions [n - 1] - positions [0];
			
			width *= -1;

			for (int d = 0; d < 3; ++d) {
				for (int k = 0; k < n; ++k) {
					for (int m = 0; m < n; ++m) {
						grids::grid::index (d, m, k) = recursion (d, m, k);
						exists (d, m, k) = true;
					}
				}
			}

			for (int k = 0; k < n; ++k) {
				for (int m = 0; m < n; ++m) {
					grids::grid::index (1, m, k) *= 2.0 / width;
					grids::grid::index (2, m, k) *= 2.0 / width * 2.0 / width;
				}
			}

			for (int d = 0; d < 3; ++d) {
				for (int k = 0; k < n; ++k) {
					grids::grid::index (d, 0, k) /= 2.0;
					grids::grid::index (d, n - 1, k) /= 2.0;
				}
			}
		}
	} /* vertical */
#else
	namespace vertical
	{
		grid::grid (axis *i_axis_ptr) :
		grids::grid (i_axis_ptr, 3) {
			TRACE ("Instantiating...");

			scale = sqrt (2.0 / (n - 1));
			pioN = std::acos (-1.0) / (n - 1);


			if (n > 1) {
				for (int i = 0; i <= n; ++i) {
					positions [i] = (i - excess_0) * (position_n - position_0) / (n - 1 - excess_n - excess_0) + position_0;
				}
			} else {
				positions [0] = (position_0 + position_n) / 2.0;
			}

			this->ood [0] = 1.0 / (positions [1] - positions [0]);
			this->ood2 [0] = 0.5 * this->ood [0];
			for (int i = 1; i < n - 1; ++i)
			{
				this->ood [i] = 1.0 / (positions [i + 1] - positions [i]);
				this->ood2 [i] = 1.0 / (positions [i + 1] - positions [i - 1]);
			}
			this->ood [n - 1] = 1.0 / (positions [n - 1] - positions [n - 2]);
			this->ood2 [n - 1] = 0.5 * this->ood [n - 1];

			// _calculate_matrix ();

			TRACE ("Instantiated...");
		}
		
		void grid::_calculate_matrix () {
			scale = 1.0 / std::sqrt (n);
			pioN = 2.0 * std::acos (-1.0) / n;

			width = positions [n] - positions [0];
			double pioL = 2.0 * std::acos (-1.0) / width;

			for (int k = 0; k < n; ++k) {
				for (int m = 0; m < n; m += 2) {
					grids::grid::index (0, m, k) = scale * std::cos (pioN * k * (m / 2)) + scale * std::cos (pioN * k * (n - m / 2));
					grids::grid::index (0, m + 1, k) = -scale * std::sin (pioN * k * (m / 2)) - scale * std::sin (pioN * k * (n - m / 2));
					grids::grid::index (1, m, k) = -scale * pioL * (m / 2) * (std::sin (pioN * k * (m / 2)) + std::sin (pioN * k * (n - m / 2)));
					grids::grid::index (1, m + 1, k) = -scale * pioL * (m / 2) * (std::cos (pioN * k * (m / 2)) + std::cos (pioN * k * (n - m / 2)));
					grids::grid::index (2, m, k) = -scale * pioL * (m / 2) * pioL * (m / 2) * (std::cos (pioN * k * (m / 2)) + std::cos (pioN * k * (n - m / 2)));
					grids::grid::index (2, m + 1, k) = scale * pioL * (m / 2) * pioL * (m / 2) * (std::sin (pioN * k * (m / 2) + std::sin (pioN * k * (n - m / 2)));
				}
			}
		}
	} /* cosine */
#endif

	namespace horizontal
	{
		grid::grid (int i_n, double i_position_0, double i_position_n, int i_excess_0, int i_excess_n, int i_ld) : 
		grids::grid (3, i_n, i_position_0, i_position_n, i_excess_0, i_excess_n, i_ld == 0 ? 2 * (i_n / 2 + 1) : i_ld) {
			TRACE ("Instantiating...");

			scale = 2.0 / sqrt (n);
			pioN = -2.0 * std::acos (-1.0) / n;
	
			if (n > 1) {
				for (int i = 0; i < n; ++i) {
					positions [i] = (i - excess_0) * (position_n - position_0) / (n - 1 - excess_n - excess_0) + position_0;
				}
			} else {
				positions [0] = (position_0 + position_n) / 2.0;
			}

			this->ood [0] = 1.0 / (positions [1] - positions [0]);
			this->ood2 [0] = 0.5 * this->ood [0];
			for (int i = 1; i < n - 1; ++i)
			{
				this->ood [i] = 1.0 / (positions [i + 1] - positions [i]);
				this->ood2 [i] = 1.0 / (positions [i + 1] - positions [i - 1]);
			}
			this->ood [n - 1] = 1.0 / (positions [n - 1] - positions [n - 2]);
			this->ood2 [n - 1] = 0.5 * this->ood [n - 1];
		
			// _calculate_matrix ();
	
			TRACE ("Instantiated...");
		}
		
		grid::grid (axis *i_axis_ptr) : 
		grids::grid (i_axis_ptr, 3, 2 * (i_axis_ptr->get_n () / 2 + 1)) {
			TRACE ("Instantiating...");

			scale = 2.0 / sqrt (n);
			pioN = -2.0 * std::acos (-1.0) / n;
	
			if (n > 1) {
				for (int i = 0; i < n; ++i) {
					positions [i] = (i - excess_0) * (position_n - position_0) / (n - 1 - excess_n - excess_0) + position_0;
				}
			} else {
				positions [0] = (position_0 + position_n) / 2.0;
			}
		
			this->ood [0] = 1.0 / (positions [1] - positions [0]);
			this->ood2 [0] = 0.5 * this->ood [0];
			for (int i = 1; i < n - 1; ++i)
			{
				this->ood [i] = 1.0 / (positions [i + 1] - positions [i]);
				this->ood2 [i] = 1.0 / (positions [i + 1] - positions [i - 1]);
			}
			this->ood [n - 1] = 1.0 / (positions [n - 1] - positions [n - 2]);
			this->ood2 [n - 1] = 0.5 * this->ood [n - 1];
		
			// _calculate_matrix ();
	
			TRACE ("Instantiated...");
		}
		
		void grid::_calculate_matrix () {
			width = positions [n - 1] - positions [0];

			for (int k = 0; k < n; ++k) {
				grids::grid::index (0, 0, k) = scale / 2.0;
				grids::grid::index (0, 1, k) = 0.0;
				grids::grid::index (1, 0, k) = 0.0;
				grids::grid::index (1, 1, k) = 0.0;
				grids::grid::index (2, 0, k) = 0.0;
				grids::grid::index (2, 1, k) = 0.0;
				for (int m = 2; m < ld; m += 2) {
					grids::grid::index (0, m, k) = scale * std::cos (pioN * k * (m / 2));
					grids::grid::index (0, m + 1, k) = scale * std::sin (pioN * k * (m / 2));
					grids::grid::index (1, m, k) = (((double) -(m / 2)) / width) * scale * std::sin (pioN * k * (m / 2));
					grids::grid::index (1, m + 1, k) = (((double) (m / 2))/ width) * scale * std::cos (pioN * k * (m / 2));
					grids::grid::index (2, m, k) = -(((double) (m / 2) * (m / 2)) / width / width) * scale * std::cos (pioN * k * (m / 2));
					grids::grid::index (2, m + 1, k) = -(((double) (m / 2) * (m / 2)) / width / width) * scale * std::sin (pioN * k * (m / 2));
				}
			}
		
			grids::grid::index (0, 1, n) = 1.0;
			for (int k = n + 1; k < ld; ++k) {
				grids::grid::index (0, k, k) = 1.0;
			}
			// for (d = 0; d < 3; ++d) {
			// 	for (k = 0; k < n; ++k) {
			// 		grids::grid::index (d, 0, k) /= 2.0;
			// 		grids::grid::index (d, n - 1, k) /= 2.0;
			// 	}
			// }
		}
		
	} /* horizontal */
} /* grids */

#endif /* end of include guard: COLLOCATION_CPP_HV4P0UOP */
