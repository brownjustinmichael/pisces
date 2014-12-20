/*!**********************************************************************
 * \file laplace.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef LAPLACE_CPP_3231A497
#define LAPLACE_CPP_3231A497

#include "laplace.hpp"

#include "linalg/exceptions.hpp"
#include "linalg-block/tridiagonal.hpp"

namespace plans
{
	template <class datatype>
	laplace_solver <datatype>::laplace_solver (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, mpi::messenger* i_messenger_ptr, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
	plans::solver <datatype> (i_element_flags, i_component_flags),
	n (i_grid_n.get_n ()), 
	ldn (i_grid_n.get_ld ()), 
	m (i_grid_m.get_n ()),
	data (i_data),
	grid_n (i_grid_n),
	grid_m (i_grid_m),
	pos_n (&i_grid_n [0]),
	pos_m (&i_grid_m [0]),
	excess_0 (grid_m.get_excess_0 ()), 
	excess_n (grid_m.get_excess_n ()),
	messenger_ptr (i_messenger_ptr) {
		TRACE ("Building laplace solver...");
		sup.resize (m * ldn);
		sub.resize (m * ldn);
		diag.resize (m * ldn);
		rhs_ptr = i_rhs;
		
		sup_ptr = &sup [0];
		sub_ptr = &sub [0];
		diag_ptr = &diag [0];
		
		supsup.resize (ldn * m);
		ipiv.resize (ldn * m);
		
		if (messenger_ptr->get_id () == 0) {
			x.resize ((m + 4 * messenger_ptr->get_np ()) * 2 * ldn);
			xipiv.resize (2 * messenger_ptr->get_np () * ldn);
		} else {
			x.resize ((m + 4) * 2 * ldn);
		}
		
		id = messenger_ptr->get_id ();
		np = messenger_ptr->get_np ();
		
		if (id != 0) {
			messenger_ptr->send (1, &pos_m [excess_0 + 1], id - 1, 0);
		}
		if (id != np - 1) {
			messenger_ptr->recv (1, &ex_pos_m, id + 1, 0);
			messenger_ptr->send (1, &pos_m [m - 1 - excess_n], id + 1, 1);
		}
		if (id != 0) {
			messenger_ptr->recv (1, &ex_pos_0, id - 1, 1);
		}
	}
	
	template <class datatype>
	void laplace_solver <datatype>::factorize () {
		TRACE ("Factorizing laplace solver...");

		double scalar = 4.0 * std::acos (-1.0) / (pos_n [n - 1] - pos_n [0]);
		int mm = m;
		int nbegin = excess_0;
		if (id != 0) {
			mm -= excess_0 + 2;
			nbegin += 1;
		}
		if (id != np - 1) {
			mm -= excess_n + 1;
		}
// #pragma omp parallel for
		for (int i = 0; i < ldn; ++i) {
			if (id != 0) {
				sub_ptr [i * m + nbegin] = 2.0 / (pos_m [nbegin + 1] - pos_m [nbegin]) / (pos_m [nbegin + 1] - ex_pos_0);
				sup_ptr [i * m + nbegin] = 2.0 / (pos_m [nbegin] - ex_pos_0) / (pos_m [nbegin + 1] - ex_pos_0);
				diag_ptr [i * m + nbegin] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [nbegin + 1] - ex_pos_0) * (1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]) + 1.0 / (pos_m [nbegin] - ex_pos_0));
			} else {
				sup_ptr [i * m + nbegin] = -1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]);
				diag_ptr [i * m + nbegin] = 1.0 / (pos_m [nbegin + 1] - pos_m [nbegin]);
				sub_ptr [i * m + nbegin] = 0.0;
			}
			for (int j = nbegin + 1; j < m - 1 - excess_n; ++j) {
				sub_ptr [i * m + j] = 2.0 / (pos_m [j + 1] - pos_m [j]) / (pos_m [j + 1] - pos_m [j - 1]);
				sup_ptr [i * m + j] = 2.0 / (pos_m [j] - pos_m [j - 1]) / (pos_m [j + 1] - pos_m [j - 1]);
				diag_ptr [i * m + j] = -scalar * (i / 2) * (i / 2) - 2.0 / (pos_m [j + 1] - pos_m [j - 1]) * (1.0 / (pos_m [j + 1] - pos_m [j]) + 1.0 / (pos_m [j] - pos_m [j - 1]));
			}
			if (id != np - 1) {
				sub_ptr [(i + 1) * m - 1 - excess_n] = 2.0 / (ex_pos_m - pos_m [m - 1 - excess_n]) / (ex_pos_m - pos_m [m - 2 - excess_n]);
				sup_ptr [(i + 1) * m - 1 - excess_n] = 2.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]) / (ex_pos_m - pos_m [m - 2 - excess_n]);
				diag_ptr [(i + 1) * m - 1 - excess_n] = -scalar * (i / 2) * (i / 2) - 2.0 / (ex_pos_m - pos_m [m - 2 - excess_n]) * (1.0 / (ex_pos_m - pos_m [m - 1 - excess_n]) + 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]));
			} else {
				sup_ptr [(i + 1) * m - 1 - excess_n] = 0.0;
				diag_ptr [(i + 1) * m - 1 - excess_n] = 1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
				sub_ptr [(i + 1) * m - 1 - excess_n] = -1.0 / (pos_m [m - 1 - excess_n] - pos_m [m - 2 - excess_n]);
			}
		}
		if (id == np - 1) {
			sup_ptr [m - 1] = 0.0;
			sub_ptr [m - 1] = 0.0;
			diag_ptr [m - 1] = 1.0;
			sup_ptr [2 * m - 1] = 0.0;
			sub_ptr [2 * m - 1] = 0.0;
			diag_ptr [2 * m - 1] = 1.0;
		}

		int info;
		
		linalg::p_block_tridiag_factorize (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], &x [0], &xipiv [0], &info, ldn, m);
	}
	
	template <class datatype>
	void laplace_solver <datatype>::execute () {
		TRACE ("Solving...");
		int mm = m;
		int nbegin = excess_0;
		if (id != 0) {
			mm -= excess_0 + 2;
			nbegin += 1;
		}
		if (id != np - 1) {
			mm -= excess_n + 1;
		}
		
		if (rhs_ptr) {
			// DEBUG (data+nbegin);
			linalg::matrix_copy (mm, ldn, rhs_ptr, data + nbegin);
		}
		
		if (id == 0) {
			linalg::scale (ldn, 0.0, data + nbegin, m);
		}
		if (id == np - 1) {
			linalg::scale (ldn, 0.0, data + m - 1 - excess_n, m);
		}
		linalg::scale (2 * m, 0.0, data);
		
		int info;
		
		// DEBUG (&sub [nbegin]);
		// DEBUG (&diag [nbegin]);
		// DEBUG (&sup [nbegin]);
		// DEBUG (&supsup [nbegin]);
		// DEBUG (&ipiv [nbegin]);
		// DEBUG (&x [0]);
		// DEBUG (&xipiv [0]);
		// DEBUG (" " << id << " " << np << " " << mm << " " << nbegin);
		linalg::p_block_tridiag_solve (id, np, mm, &sub [nbegin], &diag [nbegin], &sup [nbegin], &supsup [nbegin], &ipiv [nbegin], data + nbegin, &x [0], &xipiv [0], &info, ldn, m, m);

		// DEBUG ("DONE");
		for (int i = 0; i < ldn; ++i) {
			for (int j = nbegin - 1; j >= 0; --j) {
				data [i * m + j] = (data [i * m + j + 2] - data [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data [i * m + j + 1];
			}
			for (int j = m - excess_n; j < m; ++j) {
				data [i * m + j] = (data [i * m + j - 2] - data [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data [i * m + j - 1];
			}
		}

#ifdef CHECKNAN
		for (int j = 0; j < m; ++j) {
			for (int i = 0; i < ldn; ++i) {
				if (std::isnan (data [i * m + j])) {
					FATAL ("Nan in laplace solver.");
					throw linalg::exceptions::nan ();
				}
			}
		}
#endif

		TRACE ("Solved.");
	}
	
	template class laplace_solver <double>;
} /* plans */

#endif /* end of include guard: LAPLACE_CPP_3231A497 */
