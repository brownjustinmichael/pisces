/*!**********************************************************************
 * \file fourier.cpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2014-10-06.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef FOURIER_CPP_0B813976
#define FOURIER_CPP_0B813976

#include "linalg/utils.hpp"
#include "linalg/linalg.hpp"

#include "fourier.hpp"

namespace plans
{
	template <class datatype>
	fourier_solver <datatype>::fourier_solver (plans::grid <datatype> &i_grid_n, plans::grid <datatype> &i_grid_m, datatype& i_timestep, std::shared_ptr <plans::boundary <datatype>> i_boundary_0, std::shared_ptr <plans::boundary <datatype>> i_boundary_n, datatype *i_rhs, datatype* i_data, int *i_element_flags, int *i_component_flags) : 
	plans::solver <datatype> (i_element_flags, i_component_flags), n (i_grid_n.get_n ()), ldn (i_grid_n.get_ld ()), m (i_grid_m.get_n ()), data (i_data), timestep (i_timestep), excess_0 (i_grid_m.get_excess_0 ()), excess_n (i_grid_m.get_excess_n ()), pos_m (&i_grid_m [0]) {
		TRACE ("Building solver...");
		horizontal_matrix.resize (ldn);
		factorized_horizontal_matrix.resize (ldn);
		rhs_ptr = i_rhs;
		
		boundary_0 = i_boundary_0;
		boundary_n = i_boundary_n;
		
		ex_overlap_0 = boundary_0 ? boundary_0->get_ex_overlap () : 0;
		overlap_0 = boundary_0 ? boundary_0->get_overlap () : 0;
		ex_overlap_n = boundary_n ? boundary_n->get_ex_overlap () : 0;
		overlap_n = boundary_n ? boundary_n->get_overlap () : 0;
		lda = m + ex_overlap_n + ex_overlap_0;
		inner_m = lda - overlap_0 - overlap_n;
		
		data_temp.resize (lda * ldn);
		TRACE ("Solver built.");
	}
	
	template <class datatype>
	void fourier_solver <datatype>::factorize () {
		
		TRACE ("Factorizing...");
		
		for (int i = 0; i < ldn; ++i) {
			factorized_horizontal_matrix [i] = 1.0 + timestep * horizontal_matrix [i];
		}
		
		TRACE ("Done.");
	}
	
	template <class datatype>
	void fourier_solver <datatype>::execute () {
		TRACE ("Executing solve...");
		
		/*
			TODO Add timestep check here?
		*/
		
		// DEBUG ("Solving in n direction");
		
		linalg::scale ((ldn) * lda, 0.0, &data_temp [0]);
		
		if (rhs_ptr) {
			linalg::matrix_add_scaled (m, ldn, timestep, rhs_ptr, &data_temp [ex_overlap_0], m, lda);
		}
		
		if (boundary_0) {
			boundary_0->calculate_rhs (data + excess_0, data, &data_temp [0], &data_temp [ex_overlap_0 + excess_0], m, lda, x_solve);
		}
		if (boundary_n) {
			boundary_n->calculate_rhs (data + m - 1 - excess_n, data, &data_temp [0], &data_temp [lda - 1 - excess_n - ex_overlap_n], m, lda, x_solve);
		}
		
		linalg::matrix_add_scaled (m - 2 - excess_0 - excess_n, ldn, 1.0, data + 1 + excess_0, &data_temp [ex_overlap_0 + 1 + excess_0], m, lda);
		
		if (boundary_0) {
			boundary_0->send (&data_temp [0], lda);
		}
		if (boundary_n) {
			boundary_n->receive (&data_temp [lda - overlap_n], lda);
			boundary_n->send (&data_temp [lda - ex_overlap_n], lda);
		}
		if (boundary_0) {
			boundary_0->receive (&data_temp [ex_overlap_0], lda);
		}
		
		// DEBUG ("CHOICE RHS " << rhs_ptr [12 * m + 24] << " " << data [12 * m + 24]);
		//
		// for (int j = 0; j < lda; ++j) {
		// 	for (int i = 0; i < ldn; ++i) {
		// 		debug << data_temp [i * lda + j] << " ";
		// 	}
		// 	DEBUG ("RHS " << debug.str ());
		// 	debug.str ("");
		// }
		//
		// for (int j = 0; j < m; ++j) {
		// 	for (int i = 0; i < ldn; ++i) {
		// 		debug << data [i * m + j] << " ";
		// 	}
		// 	DEBUG ("DATA " << debug.str ());
		// 	debug.str ("");
		// }
		
		for (int j = excess_0; j < m - excess_n; ++j) {
			linalg::diagonal_solve (ldn, &factorized_horizontal_matrix [0], &data_temp [ex_overlap_0 + j], 1, lda);
		}
		linalg::matrix_copy (m, ldn, &data_temp [ex_overlap_0], data, lda);
		
		for (int i = 0; i < ldn; ++i) {
			for (int j = excess_0 - 1; j >= 0; --j) {
				data [i * m + j] = (data [i * m + j + 2] - data [i * m + j + 1]) / (pos_m [j + 2] - pos_m [j + 1]) * (pos_m [j] - pos_m [j + 1]) + data [i * m + j + 1];
			}
			for (int j = m - excess_n; j < m; ++j) {
				data [i * m + j] = (data [i * m + j - 2] - data [i * m + j - 1]) / (pos_m [j - 2] - pos_m [j - 1]) * (pos_m [j] - pos_m [j - 1]) + data [i * m + j - 1];
			}
		}
		
		TRACE ("Execution complete.");
	}
	
	template class fourier_solver <double>;
} /* plans */

#endif /* end of include guard: FOURIER_CPP_0B813976 */
