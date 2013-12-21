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
	
	enum solve_flags {
		x_solve = 0x20,
		z_solve = 0x80
	};
	
	enum initialize_flags {
		uniform_n = 0x01,
		uniform_m = 0x02
	};
	
	template <class datatype>
	class element : public bases::element <datatype>
	{
	public:
		element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags) : 
		bases::element <datatype> (i_name, 2, i_params, i_messenger_ptr, i_flags),
		axis_n (i_axis_n),
		axis_m (i_axis_m),
		n (axis_n->n), m (axis_m->n) {
			TRACE ("Instantiating...");
			
			positions [edge_n0] = i_axis_n->position_0;
			positions [edge_nn] = i_axis_n->position_n;
			positions [edge_m0] = i_axis_m->position_0;
			positions [edge_mm] = i_axis_m->position_n;
			
			excesses [edge_n0] = i_axis_n->excess_0;
			excesses [edge_nn] = i_axis_n->excess_n;
			excesses [edge_m0] = i_axis_m->excess_0;
			excesses [edge_mm] = i_axis_m->excess_n;
			
			edge_index [edge_n0] = 0;
			edge_next [edge_n0] = n;
			edge_size [edge_n0] = m;
			
			edge_index [edge_nn] = n - 1;
			edge_next [edge_nn] = n;
			edge_size [edge_nn] = m;
			
			edge_index [edge_m0] = 0;
			edge_next [edge_m0] = 1;
			edge_size [edge_m0] = n;
			
			edge_index [edge_mm] = (m - 1) * n;
			edge_next [edge_mm] = 1;
			edge_size [edge_mm] = n;
			
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
		virtual void initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
			TRACE ("Initializing " << name << "...");

			if (name == x_position) {
				initial_conditions = &(grids [0]->position ());
				flags |= uniform_m;
			} else if (name == z_position) {
				initial_conditions = &(grids [1]->position ());
				flags |= uniform_n;
			}
			// Size allowing for real FFT buffer
			scalars [name].resize (2 * (n / 2 + 1) * m, 0.0);
			if (initial_conditions) {
				if ((flags & uniform_m) && (flags & uniform_n)) {
					for (int i = 0; i < n; ++i) {
						for (int j = 0; j < m; ++j) {
							(*this) (name, i, j) = *initial_conditions;
						}
					}
				} else if (flags & uniform_m) {
					for (int j = 0; j < m; ++j) {
						utils::copy (n, initial_conditions, ptr (name, 0, j), 1, m);
					}
				} else if (flags & uniform_n) {
					for (int i = 0; i < n; ++i) {
						utils::copy (m, initial_conditions, ptr (name, i, 0));
					}
				} else {
					utils::copy (n * m, initial_conditions, ptr (name));
				}
			}
			
			TRACE ("Initialized.");
		}
		
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::explicit_reset ()
		 *********************************************************************/
		inline void explicit_reset () {
			bases::element <datatype>::explicit_reset ();
			
			TRACE ("Explicit reset start...");
			
			for (typename std::map <int, std::vector <datatype> >::iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				if (iter->first < 0) {
					utils::scale ((2 * (n / 2 + 1)) * m, 0.0, &(iter->second) [0]);
				}
			}
			if (flags & z_solve) {
				flags &= ~z_solve;
				flags |= x_solve;
			} else {
				flags &= ~x_solve;
				flags |= z_solve;
			}
			
			TRACE ("Explicit reset end.");
		}
	
	protected:
		using bases::element <datatype>::scalars;
		using bases::element <datatype>::grids;
		using bases::element <datatype>::name;
		using bases::element <datatype>::failsafe_dump;
		using bases::element <datatype>::messenger_ptr;
		using bases::element <datatype>::flags;
		
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
				element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags) : 
				two_d::element <datatype> (i_axis_n, i_axis_m, i_name, i_params, i_messenger_ptr, i_flags) {
					TRACE ("Instantiating...");
					
					two_d::element <datatype>::set_grid (new bases::fourier::grid <datatype> (axis_n));
					two_d::element <datatype>::set_grid (new bases::chebyshev::grid <datatype> (axis_m), 1);
					initialize (x_position);
					initialize (z_position);
					
					TRACE ("Instantiated.");
				}
				virtual ~element () {}
				
				virtual void initialize (int name, datatype* initial_conditions = NULL, int flags = 0x00) {
					TRACE ("Initializing...");
					two_d::element <datatype>::initialize (name, initial_conditions, flags);
					if ((name != x_position) && (name != z_position)) {
					    if (flags & only_forward_horizontal) {
							element <datatype>::add_forward_horizontal_transform (new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name)));
						} else if (!(flags & no_transform)) {
							DEBUG ("Adding all transforms for " << name);
			   				element <datatype>::add_forward_horizontal_transform (new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name)));
							element <datatype>::add_inverse_horizontal_transform (new horizontal_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse));
							element <datatype>::add_forward_vertical_transform (new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name)));
							element <datatype>::add_inverse_vertical_transform (new vertical_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, inverse));
						}
					}
				}
		
			protected:
				using two_d::element <datatype>::initialize;
				using two_d::element <datatype>::positions;
				using two_d::element <datatype>::excesses;
				using two_d::element <datatype>::n;
				using two_d::element <datatype>::m;
				using two_d::element <datatype>::axis_n;
				using two_d::element <datatype>::axis_m;
				using two_d::element <datatype>::params;
				using two_d::element <datatype>::ptr;
				using two_d::element <datatype>::grids;
			};
			
			template <class datatype>
			class advection_diffusion_element : public element <datatype>
			{
			public:
				advection_diffusion_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags);
				
				virtual ~advection_diffusion_element () {}
			
				datatype calculate_timestep ();
			
			private:
				using element <datatype>::flags;
				using element <datatype>::params;
				using element <datatype>::initialize;
				using element <datatype>::n;
				using element <datatype>::m;
				using element <datatype>::name;
				using element <datatype>::normal_stream;
				using element <datatype>::transform_stream;
				using element <datatype>::cell_n;
				using element <datatype>::cell_m;
				using element <datatype>::grids;
				using element <datatype>::ptr;
				using element <datatype>::matrix_ptr;
				using element <datatype>::messenger_ptr;
				using element <datatype>::timestep;
				using element <datatype>::alpha_0;
				using element <datatype>::alpha_n;
				using element <datatype>::solvers;
			};
			
			template <class datatype>
			class convection_element : public element <datatype>
			{
			public:
				convection_element (bases::axis *i_axis_n, bases::axis *i_axis_m, int i_name, io::parameters <datatype>& i_params, bases::messenger* i_messenger_ptr, int i_flags);
				
				virtual ~convection_element () {}
			
				datatype calculate_timestep ();
			
			private:
				using element <datatype>::flags;
				using element <datatype>::params;
				using element <datatype>::initialize;
				using element <datatype>::n;
				using element <datatype>::m;
				using element <datatype>::name;
				using element <datatype>::normal_stream;
				using element <datatype>::transform_stream;
				using element <datatype>::cell_n;
				using element <datatype>::cell_m;
				using element <datatype>::grids;
				using element <datatype>::ptr;
				using element <datatype>::matrix_ptr;
				using element <datatype>::messenger_ptr;
				using element <datatype>::timestep;
				using element <datatype>::alpha_0;
				using element <datatype>::alpha_n;
				using element <datatype>::solvers;
			};
		} /* chebyshev */
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
