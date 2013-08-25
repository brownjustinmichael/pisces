/*!***********************************************************************
 * \file collocation_one_d.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cassert>
#include <cmath>
#include <iomanip>
#include <exception>
#include "../config.hpp"
#include "../bases/collocation.hpp"
#include "collocation_one_d.hpp"

namespace one_d
{
	namespace chebyshev
	{
		template <class datatype>
		grid <datatype>::grid (int i_M, int i_N, datatype i_scale, datatype i_width) : 
		bases::collocation_grid <datatype> (3, i_M, i_N) {
			int d, m, k;
			scale = i_scale;
			width = i_width;
			pioN = std::acos (-1.0) / (i_N - 1);
			exists_array.resize (i_M * i_N * 3, false);
	
			TRACE ("Instantiating...");
	
			for (d = 0; d < 3; ++d) {
				for (k = 0; k < i_N; ++k) {
					for (m = 0; m < i_M; ++m) {
						bases::collocation_grid <datatype>::index (d, m, k) = recursion (d, m, k);
						exists (d, m, k) = true;
					}
				}
			}
	
			for (k = 0; k < i_N; ++k) {
				for (m = 0; m < i_M; ++m) {
					bases::collocation_grid <datatype>::index (1, m, k) *= 2.0 / width;
					bases::collocation_grid <datatype>::index (2, m, k) *= 2.0 / width * 2.0 / width;
				}
			}
	
			for (d = 0; d < 3; ++d) {
				for (k = 0; k < i_N; ++k) {
					bases::collocation_grid <datatype>::index (d, 0, k) /= 2.0;
					bases::collocation_grid <datatype>::index (d, i_M - 1, k) /= 2.0;
				}
			}
	
			TRACE ("Instantiated...");
		}

		template <class datatype>
		datatype grid <datatype>::recursion (int d, int m, int k) {		
			assert (d >= 0);
			assert (d < 3);
			assert (m >= 0);
			assert (k >= 0);

			// Use the recursion relations to calculate the correct value of the polynomial at the collocation point
			if (exists (d, m, k)) {
				// This value has already been calculated
				return bases::collocation_grid <datatype>::index (d, m, k);
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
		
	template class grid <double>;
	template class grid <float>;
	} /* chebyshev */
	
	namespace fourier
	{
		template <class datatype>
		grid <datatype>::grid (int i_M, int i_N, datatype i_scale, datatype i_width) : 
		bases::collocation_grid <datatype> (3, i_M, i_N) {
			int d, m, k;
			scale = i_scale;
			width = i_width;
			pioN = std::acos (-1.0) / (i_N - 1);
	
			TRACE ("Instantiating...");
	
			for (k = 0; k < i_N; ++k) {
				for (m = 0; m < i_M; ++m) {
					bases::collocation_grid <datatype>::index (0, m, k) = scale * std::cos (pioN * k * m);
					bases::collocation_grid <datatype>::index (1, m, k) = (((datatype) -m) * 2.0 / width) * scale * std::sin (pioN * k * m);
					bases::collocation_grid <datatype>::index (2, m, k) = -(((datatype) m * m) * 4.0 / width / width) * scale * std::cos (pioN * k * m);
				}
			}
	
			for (d = 0; d < 3; ++d) {
				for (k = 0; k < i_N; ++k) {
					bases::collocation_grid <datatype>::index (d, 0, k) /= 2.0;
					bases::collocation_grid <datatype>::index (d, i_M - 1, k) /= 2.0;
				}
			}
	
			TRACE ("Instantiated...");
		}
		
	template class grid <double>;
	template class grid <float>;
	} /* fourier */
} /* one_d */
	