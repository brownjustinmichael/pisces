/*!**********************************************************************
 * \file element_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_TWO_D_HPP_CJ68F4IB
#define ELEMENT_TWO_D_HPP_CJ68F4IB


#include "../bases/element.hpp"
#include <cmath>
#include <sstream>
#include "transform_two_d.hpp"
#include "../utils/io.hpp"
#include "../config.hpp"

namespace two_d
{
	enum edges {
		edge_n0 = 0, // Start at 0, 0, increment by n
		edge_nn = 1, // Start at n, 0, increment by n
		edge_m0 = 2, // Start at 0, 0, increment by 1
		edge_mm = 3 // Start at 0, m, increment by 1
	};
	
	enum solve_element_flags {
		x_solve = 0x20,
		z_solve = 0x80
	};
	
	enum initialize_element_flags {
		uniform_n = 0x01,
		uniform_m = 0x02
	};
	
	template <class datatype>
	class element : public bases::element <datatype>
	{
	public:
		element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags) : 
		bases::element <datatype> (i_name, 2, i_params, i_messenger_ptr, i_element_flags),
		axis_n (i_axis_n),
		axis_m (i_axis_m),
		n (axis_n->n), m (axis_m->n) {
			TRACE ("Instantiating...");
			
			if (i_messenger_ptr->get_id () == 0) {
				alpha_0 = 0.0;
			} else {
				alpha_0 = 0.5;
			}
			if (i_messenger_ptr->get_id () == i_messenger_ptr->get_np () - 1) {
				alpha_n = 0.0;
			} else {
				alpha_n = 0.5;
			}
			
			cell_n.resize (n * m);
			cell_m.resize (n * m);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					cell_n [i * m + j] = i;
					cell_m [i * m + j] = j;
				}
			}
			
			std::ostringstream convert;
			convert << name;
			failsafe_dump.reset (new io::output (new io::two_d::netcdf (n, m), "dump.dat"));
			failsafe_dump->template append <int> ("i", &cell_n [0]);
			failsafe_dump->template append <int> ("j", &cell_m [0]);
			
			max_timestep = i_params.get <datatype> ("time.max");
			timestep_safety = i_params.get <datatype> ("time.safety");
			
			TRACE ("Instantiated.");
		}
		
		virtual ~element () {}
		
		/*!*******************************************************************
		 * \brief Get the datatype reference to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * \return A datatype reference to the first element of the named scalar
		 *********************************************************************/
		inline datatype& operator[] (int name) {
			if (scalars.find (name) == scalars.end ()) {
				FATAL ("Index " << name << " not found in element.");
				throw 0;
			}
			return scalars [name] [0];
		}
		
		inline datatype& operator() (int name, int i = 0, int j = 0) {
			return bases::element <datatype>::operator() (name, i * m + j);
		}
		
		inline datatype* ptr (int name, int i = 0) {
			return bases::element <datatype>::ptr (name, i);
		}

		inline datatype* ptr (int name, int i, int j) {
			return bases::element <datatype>::ptr (name, i * m + j);
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");

			if (name == x_position) {
				initial_conditions = &((*grids [0]) [0]);
				i_flags |= uniform_m;
			} else if (name == z_position) {
				initial_conditions = &((*grids [1]) [0]);
				i_flags |= uniform_n;
			}
			// Size allowing for real FFT buffer
			scalars [name].resize (grids [0]->ld * m, 0.0);
			if (initial_conditions) {
				if ((i_flags & uniform_m) && (i_flags & uniform_n)) {
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							(*this) (name, i, j) = *initial_conditions;
						}
					}
				} else if (i_flags & uniform_m) {
					for (int j = 0; j < m; ++j) {
						utils::copy (n, initial_conditions, ptr (name, 0, j), 1, m);
					}
				} else if (i_flags & uniform_n) {
					for (int i = 0; i < n; ++i) {
						utils::copy (m, initial_conditions, ptr (name, i, 0));
					}
				} else {
					utils::copy (n * m, initial_conditions, ptr (name));
				}
			}
			
			TRACE ("Initialized.");
			return this->ptr (name);
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset ();
			
			TRACE ("Explicit reset start...");
			
			if (element_flags [state] & z_solve) {
				element_flags [state] &= ~z_solve;
				element_flags [state] |= x_solve;
			} else {
				element_flags [state] &= ~x_solve;
				element_flags [state] |= z_solve;
			}
			
			/*
				TODO It'd be nice if this could be moved to a different method
			*/
			
			TRACE ("Explicit reset end.");
		}
		
		virtual datatype calculate_timestep (int i, int j, datatype *i_position = NULL, datatype *i_velocity = NULL, int flags = 0x00) = 0;
		
		inline datatype calculate_min_timestep (int i_m = -1, datatype *i_position = NULL, datatype *i_velocity = NULL, int flags = 0x00) {
			double shared_min = max_timestep;
			#pragma omp parallel 
			{
				double min = std::numeric_limits <double>::max ();
				int nmin = 1, nmax = n - 1;
				if (flags & profile_timestep) {
					nmin = 0;
					nmax = 1;
				}
				if (i_m == -1) {
					i_m = m;
				}
				#pragma omp for nowait
					for (int j = 1; j < i_m - 1; ++j) {
						for (int i = nmin; i < nmax; ++i) {
							min = std::min (calculate_timestep (i, j, i_position, i_velocity, flags), min);
						}
					}
				#pragma omp critical 
				{
					shared_min = std::min (shared_min, min);
				}
			}
			if (shared_min < timestep || shared_min > 2.0 * timestep) {
				return shared_min * timestep_safety;
			} else {
				return timestep;
			}
		}
	
	protected:
		using bases::element <datatype>::scalars;
		using bases::element <datatype>::grids;
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		using bases::element <datatype>::element_flags;
		using bases::element <datatype>::timestep;
		
		bases::axis *axis_n, *axis_m;
		int &n; //!< The number of elements in each 1D array
		int &m;
		datatype alpha_0, alpha_n;
		std::map <int, datatype> positions; //!< A vector of the edge positions
		std::map <int, int> excesses; //!< A vector of the edge positions
		std::vector<int> cell_n; //!< An integer array for tracking each cell number for output
		std::vector<int> cell_m; //!< An integer array for tracking each cell number for output
		
		std::map <int, int> edge_index;
		std::map <int, int> edge_next;
		std::map <int, int> edge_size;
		
		datatype max_timestep;
		datatype timestep_safety;
	};
	
	namespace fourier
	{
		namespace chebyshev
		{
			template <class datatype>
			class element : public two_d::element <datatype>
			{
			public:
				/*!*******************************************************************
				 * \copydoc one_d::element::element ()
				 *********************************************************************/
				element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags) : 
				two_d::element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_element_flags) {
					TRACE ("Instantiating...");
					initialize (x_position, "x");
					initialize (z_position, "z");
					if (i_params ["parallel.transform.subthreads"].IsDefined ()) {
						transform_threads = i_params.get <int> ("parallel.transform.subthreads");
					} else {
						transform_threads = 0;
					}
					TRACE ("Instantiated.");
				}
				virtual ~element () {}
				
				int &get_mode () {
					return mode;
				}
				
				virtual std::shared_ptr <bases::grid <datatype>> generate_grid (int index = 0) {
					if (index == 0) {
						return std::shared_ptr <bases::grid <datatype>> (new typename horizontal_grid <datatype>::type (axis_n));
					} else if (index == 1) {
						return std::shared_ptr <bases::grid <datatype>> (new typename vertical_grid <datatype>::type (axis_m));
					} else {
						throw 0;
					}
				}
				
				virtual datatype *_initialize (int name, int i_threads = 0, datatype* initial_conditions = NULL, int i_element_flags = 0x00) {
					TRACE ("Initializing...");
					two_d::element <datatype>::_initialize (name, initial_conditions, i_element_flags);
					/*
						TODO Fix flaggin
					*/
					if ((name != x_position) && (name != z_position)) {
						element <datatype>::add_transform (name, new master_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal , &element_flags [state], &element_flags [name], transform_threads));
						// 					    if (i_element_flags & only_forward_horizontal) {
						// 	element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], i_threads), forward_horizontal);
						// } else if (!(i_element_flags & no_transform)) {
						// 			   				element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], i_threads), forward_horizontal);
						// 	element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse, &element_flags [state], &element_flags [name], i_threads), inverse_horizontal);
						// 	element <datatype>::add_transform (name, new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], i_threads), forward_vertical);
						// 	element <datatype>::add_transform (name, new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse, &element_flags [state], &element_flags [name], i_threads), inverse_vertical);
						// }
					}
					return this->ptr (name);
				}
		
				using two_d::element <datatype>::ptr;
			protected:
				using two_d::element <datatype>::initialize;
				using two_d::element <datatype>::element_flags;
				using two_d::element <datatype>::positions;
				using two_d::element <datatype>::excesses;
				using two_d::element <datatype>::n;
				using two_d::element <datatype>::m;
				using two_d::element <datatype>::axis_n;
				using two_d::element <datatype>::axis_m;
				using two_d::element <datatype>::params;
				using two_d::element <datatype>::grids;
				
				int transform_threads;
				static int mode;
			};
		}
		
		namespace cosine
		{	
			template <class datatype>
			class element : public two_d::element <datatype>
			{
			public:
				/*!*******************************************************************
				 * \copydoc one_d::element::element ()
				 *********************************************************************/
				element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters& i_params, bases::messenger* i_messenger_ptr, int i_element_flags) : 
				two_d::element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_element_flags) {
					TRACE ("Instantiating...");
					initialize (x_position, "x");
					initialize (z_position, "z");
					
					if (i_params ["parallel.transform.subthreads"].IsDefined ()) {
						transform_threads = i_params.get <int> ("parallel.transform.subthreads");
					} else {
						transform_threads = 0;
					}
					TRACE ("Instantiated.");
				}
				
				virtual ~element () {}
				
				int &get_mode () {
					return mode;
				}
				
				virtual std::shared_ptr <bases::grid <datatype>> generate_grid (int index = 0) {
					if (index == 0) {
						return std::shared_ptr <bases::grid <datatype>> (new typename horizontal_grid <datatype>::type (axis_n));
					} else if (index == 1) {
						return std::shared_ptr <bases::grid <datatype>> (new typename vertical_grid <datatype>::type (axis_m));
					} else {
						throw 0;
					}
				}
				
				virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int i_element_flags = 0x00) {
					TRACE ("Initializing...");
					two_d::element <datatype>::_initialize (name, initial_conditions, i_element_flags);
					/*
						TODO Fix flaggin
					*/
					if ((name != x_position) && (name != z_position)) {
						bases::element <datatype>::add_transform (name, new master_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal, &element_flags [state], &element_flags [name], transform_threads));
						// 					    if (i_element_flags & only_forward_horizontal) {
						// 	element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], transform_threads), forward_horizontal);
						// } else if (!(i_element_flags & no_transform)) {
						// 			   				element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], transform_threads), forward_horizontal);
						// 	element <datatype>::add_transform (name, new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse, &element_flags [state], &element_flags [name], transform_threads), inverse_horizontal);
						// 	element <datatype>::add_transform (name, new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, NULL, &element_flags [state], &element_flags [name], transform_threads), forward_vertical);
						// 	element <datatype>::add_transform (name, new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse, &element_flags [state], &element_flags [name], transform_threads), inverse_vertical);
						// }
					}
					return this->ptr (name);
				}
		
				using two_d::element <datatype>::ptr;
			protected:
				using two_d::element <datatype>::initialize;
				using two_d::element <datatype>::element_flags;
				using two_d::element <datatype>::positions;
				using two_d::element <datatype>::excesses;
				using two_d::element <datatype>::n;
				using two_d::element <datatype>::m;
				using two_d::element <datatype>::axis_n;
				using two_d::element <datatype>::axis_m;
				using two_d::element <datatype>::params;
				using two_d::element <datatype>::grids;
				
				int transform_threads;
				
				static int mode;
			};
		}
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
