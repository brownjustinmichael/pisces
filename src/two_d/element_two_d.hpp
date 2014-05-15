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
#include "../utils/rezone.hpp"
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
		element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
		bases::element <datatype> (i_name, 2, i_params, i_messenger_ptr, i_element_flags),
		n (i_axis_n.get_n ()), m (i_axis_m.get_n ()) {
			TRACE ("Instantiating...");
			axes [0] = i_axis_n;
			axes [1] = i_axis_m;
			
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
		
		virtual void setup_profile (std::shared_ptr <io::output> output_stream, int flags = 0x00) {
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				output_stream->template append <datatype> (scalar_names [iter->first], ptr (iter->first));
				output_stream->template append_functor <datatype> ("rms_" + scalar_names [iter->first], new io::root_mean_square_functor <datatype> (ptr (iter->first), n, m));
				output_stream->template append_functor <datatype> ("avg_" + scalar_names [iter->first], new io::average_functor <datatype> (ptr (iter->first), n, m));
			}
			output_stream->template append_scalar <datatype> ("t", &duration);
			output_stream->template append_scalar <const int> ("mode", &(get_mode ()));
			
			bases::element <datatype>::setup_profile (output_stream, flags);
		}
	
		/*!*******************************************************************
		 * \copydoc bases::element <datatype>::initialize ()
		 *********************************************************************/
		virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			// Size allowing for real FFT buffer
			scalars [name].resize (grids [0]->get_ld () * m, 0.0);
			if (name == x_position) {
				for (int j = 0; j < m; ++j) {
					utils::copy (n, &((*grids [0]) [0]), ptr (name, 0, j), 1, m);
				}
			} else if (name == z_position) {
				for (int i = 0; i < n; ++i) {
					utils::copy (m, &((*grids [1]) [0]), ptr (name, i, 0));
				}
			} else {
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
			}
			
			TRACE ("Initialized.");
			return this->ptr (name);
		}
		
		virtual datatype calculate_timestep (int i, int j, io::virtual_dump *dump = NULL) = 0;
		
		inline datatype calculate_min_timestep (io::virtual_dump *dump = NULL) {
			double shared_min = max_timestep;
			#pragma omp parallel 
			{
				double min = std::numeric_limits <double>::max ();
				#pragma omp for nowait
					for (int j = 1; j < (dump ? dump->dims ["z"] [1] : m) - 1; ++j) {
						for (int i = 0; i < (dump ? dump->dims ["z"] [0] : n); ++i) {
							min = std::min (calculate_timestep (i, j, dump), min);
						}
					}
				#pragma omp critical 
				{
					shared_min = std::min (shared_min, min);
				}
			}
			if ((shared_min < timestep || shared_min > 2.0 * timestep) || dump) {
				return shared_min * timestep_safety;
			} else {
				return timestep;
			}
		}
		
		virtual std::shared_ptr <io::virtual_dump> make_dump (int flags = 0x00) {
			std::shared_ptr <io::virtual_dump> dump (new io::virtual_dump);
			
			std::shared_ptr <io::output> virtual_output (new io::output (new io::two_d::virtual_format (&*dump, n, m)));
			bases::element <datatype>::setup_output (virtual_output);
			
			virtual_output->to_file ();
			return dump;
		}
		
		virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) = 0;
		
		virtual std::shared_ptr <io::virtual_dump> make_rezoned_dump (datatype *positions, io::virtual_dump *old_dump, int flags = 0x00) {
			std::shared_ptr <io::virtual_dump> dump (new io::virtual_dump);
			
			bases::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <bases::grid <datatype>> vertical_grid = generate_grid (&vertical_axis);
			
			utils::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, old_dump, &*dump);
			
			return dump;
		}
		
		/*
			TODO Combine these?
		*/
		
		virtual void get_zoning_positions (datatype *positions) {
			if (messenger_ptr->get_id () == 0) {
				datatype temp [messenger_ptr->get_np () * 2];
				temp [0] = axes [1].get_position_0 ();
				temp [1] = axes [1].get_position_n ();
				messenger_ptr->template gather <datatype> (2, temp, temp);
				for (int i = 0; i < messenger_ptr->get_np (); ++i) {
					positions [i] = temp [2 * i];
				}
				positions [messenger_ptr->get_np ()] = temp [messenger_ptr->get_np () * 2 - 1];
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			} else {
				datatype temp [2] = {axes [1].get_position_0 (), axes [1].get_position_n ()};
				messenger_ptr->template gather <datatype> (2, temp, temp);
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			}
		}
	
	protected:
		using bases::element <datatype>::scalars;
		using bases::element <datatype>::grids;
		using bases::element <datatype>::name;
		using bases::element <datatype>::messenger_ptr;
		using bases::element <datatype>::element_flags;
		using bases::element <datatype>::timestep;
		using bases::element <datatype>::axes;
		using bases::element <datatype>::scalar_names;
		using bases::element <datatype>::get_mode;
		using bases::element <datatype>::duration;
		
		int n; //!< The number of elements in each 1D array
		int m;
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
				element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
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
				
				const int &get_mode () {
					return mode;
				}
				
				virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) {
					if (index == 0) {
						return std::shared_ptr <bases::grid <datatype>> (new typename horizontal_grid <datatype>::type (axis));
					} else if (index == 1 || index == -1) {
						return std::shared_ptr <bases::grid <datatype>> (new typename vertical_grid <datatype>::type (axis));
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
						element <datatype>::add_transform (name, std::shared_ptr <master_transform <datatype> > (new master_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal , &element_flags [state], &element_flags [name], transform_threads)));
					}
					TRACE ("Initialized.");
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
				element (bases::axis i_axis_n, bases::axis i_axis_m, int i_name, io::parameters& i_params, utils::messenger* i_messenger_ptr, int i_element_flags) : 
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
				
				const int &get_mode () {
					return mode;
				}
				
				virtual std::shared_ptr <bases::grid <datatype>> generate_grid (bases::axis *axis, int index = -1) {
					if (index == 0) {
						return std::shared_ptr <bases::grid <datatype>> (new typename horizontal_grid <datatype>::type (axis));
					} else if (index == 1 || index == -1) {
						return std::shared_ptr <bases::grid <datatype>> (new typename vertical_grid <datatype>::type (axis));
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
						bases::element <datatype>::add_transform (name, std::shared_ptr <master_transform <datatype> > (new master_transform <datatype> (*grids [0], *grids [1], ptr (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal, &element_flags [state], &element_flags [name], transform_threads)));
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
				using two_d::element <datatype>::params;
				using two_d::element <datatype>::grids;
				
				int transform_threads;
				
				static int mode;
			};
		}
	} /* fourier */
} /* two_d */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
