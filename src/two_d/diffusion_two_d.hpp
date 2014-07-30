/*!**********************************************************************
 * \file diffusion_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2013-10-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DIFFUSION_TWO_D_HPP_YHECX9VS
#define DIFFUSION_TWO_D_HPP_YHECX9VS

#include "../bases/plan.hpp"
#include "solver_two_d.hpp"

/*
	TODO Should non-axis diffusion be in the explicit or implicit rhs?
*/

namespace two_d
{
	namespace fourier
	{
		// namespace chebyshev
		// {
			template <class datatype>
			class horizontal_diffusion : public implicit_plan <datatype>
			{
			public:
				horizontal_diffusion (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
				implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
				coeff (i_coeff),
				alpha (i_alpha) {
					TRACE ("Instantiating...");
					pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
					if (matrix_n) {
						for (int i = 0; i < ldn; ++i) {
							matrix_n [i] = coeff * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
						}
					} else {
						WARN ("No matrix");
					}
				}
				
				horizontal_diffusion (bases::master_solver <datatype> &i_solver, datatype i_coeff, datatype i_alpha) :
				implicit_plan <datatype> (i_solver),
				coeff (i_coeff),
				alpha (i_alpha) {
					TRACE ("Instantiating...");
					pioL2 = 4.0 * (std::acos (-1.0) * std::acos (-1.0) / (grid_n [n - 1] - grid_n [0]) / (grid_n [n - 1] - grid_n [0]));
					if (matrix_n) {
						for (int i = 0; i < ldn; ++i) {
							matrix_n [i] = coeff * alpha * pioL2 * (datatype) ((i / 2) * (i / 2));
						}
					} else {
						WARN ("No matrix");
					}
				}
				
				virtual ~horizontal_diffusion () {}
			
				void execute () {	
					TRACE ("Operating..." << element_flags);
					std::stringstream debug;
					datatype max = 0.0;
					int index = 0;
					for (int j = 0; j < m; ++j) {
						for (int i = 0; i < ldn; ++i) {
							if (data_out [i * m + j] < 0.0) {
								if (data_out [i * m + j] < -max) {
									max = -data_out [i * m + j];
									index = i * m + j;
								}
							} else {
								if (data_out [i * m + j] > max) {
									max = data_out [i * m + j];
									index = i * m + j;
								}
							}
						}
					}
					DEBUG ("MAX BEFORE HDIFF " << (index / m) << " " << (index % m) << " " << data_out [index] << " " << data_in [index]);
					DEBUG ("CHOICE FROM HDIFF " << data_out [12 * m + 24] << " " << data_in [12 * m + 24]);
					if (*component_flags & x_solve) {
						if (1.0 - alpha != 0.0) {
							#pragma omp parallel for
							for (int i = 0; i < ldn; ++i) {
								utils::add_scaled (m, -coeff * (1.0 - alpha) * pioL2 * (i / 2) * (i / 2), data_in + i * m, data_out + i * m);
							}
						}
					} else {
						#pragma omp parallel for
						for (int i = 0; i < ldn; ++i) {
							utils::add_scaled (m, -coeff * pioL2 * (i / 2) * (i / 2), data_in + i * m, data_out + i * m);
						}
					}
					TRACE ("Operation complete.");
					
					max = 0.0;
					index = 0;
					for (int j = 0; j < m; ++j) {
						for (int i = 0; i < ldn; ++i) {
							if (data_out [i * m + j] < 0.0) {
								if (data_out [i * m + j] < -max) {
									max = -data_out [i * m + j];
									index = i * m + j;
								}
							} else {
								if (data_out [i * m + j] > max) {
									max = data_out [i * m + j];
									index = i * m + j;
								}
							}
						}
					}
					DEBUG ("MAX FROM HDIFF " << (index / m) << " " << (index % m) << " " << data_out [index] << " " << data_in [index]);
					DEBUG ("CHOICE FROM HDIFF " << data_out [12 * m + 24] << " " << data_in [12 * m + 24]);
				}
				
				using implicit_plan <datatype>::element_flags;
				using implicit_plan <datatype>::component_flags;
			
			private:
				datatype coeff;
				datatype alpha;
				datatype pioL2;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::ldn;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::grid_n;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
			};
			
			template <class datatype>
			class vertical_diffusion : public implicit_plan <datatype>
			{
			public:
				vertical_diffusion (bases::grid <datatype> &i_grid_n, bases::grid <datatype> &i_grid_m, datatype i_coeff, datatype i_alpha, datatype *i_matrix_n, datatype *i_matrix_m, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
				implicit_plan <datatype> (i_grid_n, i_grid_m, i_matrix_n, i_matrix_m, i_data_in, i_data_out, i_element_flags, i_component_flags),
				coeff (i_coeff),
				alpha (i_alpha) {
					if (matrix_m) {
						for (int j = 0; j < m; ++j) {
							utils::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
						}
					} else {
						WARN ("No matrix");
					}
				}
				
				vertical_diffusion (bases::master_solver <datatype> &i_solver, datatype i_coeff, datatype i_alpha) :
				implicit_plan <datatype> (i_solver),
				coeff (i_coeff),
				alpha (i_alpha) {
					if (matrix_m) {
						for (int j = 0; j < m; ++j) {
							utils::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
						}
					} else {
						WARN ("No matrix");
					}
				}
				
				virtual ~vertical_diffusion () {}
			
				void execute () {	
					TRACE ("Operating..." << element_flags);
					
					datatype max = 0.0;
					int index = 0;
					for (int j = 0; j < m; ++j) {
						for (int i = 0; i < ldn; ++i) {
							if (data_out [i * m + j] < 0.0) {
								if (data_out [i * m + j] < -max) {
									max = -data_out [i * m + j];
									index = i * m + j;
								}
							} else {
								if (data_out [i * m + j] > max) {
									max = data_out [i * m + j];
									index = i * m + j;
								}
							}
						}
					}
					DEBUG ("MAX BEFORE VDIFF " << (index / m) << " " << (index % m) << " " << data_out [index] << " " << data_in [index]);
					DEBUG ("CHOICE BEFORE VDIFF " << data_out [12 * m + 24] << " " << data_in [12 * m + 24]);
					if (*component_flags & z_solve) {
						if (1.0 - alpha != 0.0) {
							utils::matrix_matrix_multiply (m, ldn, m, coeff * (1.0 - alpha), grid_m.get_data (2), data_in, 1.0, data_out, m);
						}
					} else {
						utils::matrix_matrix_multiply (m, ldn, m, coeff, grid_m.get_data (2), data_in, 1.0, data_out, m);
					}
					
					max = 0.0;
					index = 0;
					for (int j = 0; j < m; ++j) {
						for (int i = 0; i < ldn; ++i) {
							if (data_out [i * m + j] < 0.0) {
								if (data_out [i * m + j] < -max) {
									max = -data_out [i * m + j];
									index = i * m + j;
								}
							} else {
								if (data_out [i * m + j] > max) {
									max = data_out [i * m + j];
									index = i * m + j;
								}
							}
						}
					}
					DEBUG ("MAX FROM VDIFF " << (index / m) << " " << (index % m) << " " << data_out [index] << " " << data_in [index]);
					DEBUG ("CHOICE FROM VDIFF " << data_out [12 * m + 24] << " " << data_in [12 * m + 24]);
					
					TRACE ("Operation complete.");
				}
			
				using implicit_plan <datatype>::element_flags;
				using implicit_plan <datatype>::component_flags;
			
			private:
				datatype coeff;
				datatype alpha;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::ldn;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
				using implicit_plan <datatype>::matrix_m;
				using implicit_plan <datatype>::grid_n;
				using implicit_plan <datatype>::grid_m;
			};
			
			template <class datatype>
			class finite_vertical_diffusion : public implicit_plan <datatype>
			{
			public:
				finite_vertical_diffusion (bases::master_solver <datatype> &i_solver, datatype i_coeff, datatype i_alpha) :
				implicit_plan <datatype> (i_solver),
				coeff (i_coeff),
				alpha (i_alpha) {
					TRACE ("Instantiating...");
					for (int j = 0; j < m; ++j) {
						utils::add_scaled (m, -coeff * alpha, grid_m.get_data (2) + j, matrix_m + j, m, m);
					}
				}
				
				virtual ~finite_vertical_diffusion () {}
			
				void execute () {	
					TRACE ("Operating...");
					
					if (*component_flags & z_solve) {
						if (1.0 - alpha != 0.0) {
							for (int i = 1; i < m - 1; ++i) {
								datatype pos_m1 = grid_m [i - 1];
								datatype pos_0 = grid_m [i];
								datatype pos_1 = grid_m [i + 1];
								for (int j = 0; j < ldn; ++j) {
									data_out [j * m + i] += coeff * (1.0 - alpha) * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
								}
							}
						}
					} else {
						for (int i = 1; i < m - 1; ++i) {
							datatype pos_m1 = grid_m [i - 1];
							datatype pos_0 = grid_m [i];
							datatype pos_1 = grid_m [i + 1];
							for (int j = 0; j < ldn; ++j) {
								data_out [j * m + i] += coeff * 2.0 * ((data_in [j * m + i + 1] - data_in [j * m + i]) / (pos_1 - pos_0) - (data_in [j * m + i] - data_in [j * m + i - 1]) / (pos_0 - pos_m1)) / (pos_1 - pos_m1);
							}
						}
					}
					
					TRACE ("Operation complete.");
				}
			
				using implicit_plan <datatype>::element_flags;
				using implicit_plan <datatype>::component_flags;
			
			private:
				datatype coeff;
				datatype alpha;
				using implicit_plan <datatype>::n;
				using implicit_plan <datatype>::ldn;
				using implicit_plan <datatype>::m;
				using implicit_plan <datatype>::data_in;
				using implicit_plan <datatype>::data_out;
				using implicit_plan <datatype>::matrix_n;
				using implicit_plan <datatype>::matrix_m;
				using implicit_plan <datatype>::grid_n;
				using implicit_plan <datatype>::grid_m;
			};
		// } /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: DIFFUSION_TWO_D_HPP_YHECX9VS */
