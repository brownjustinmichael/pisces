/*!**********************************************************************
 * \file element_two_d.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_TWO_D_HPP_CJ68F4IB
#define ELEMENT_TWO_D_HPP_CJ68F4IB


#include <cmath>
#include <sstream>

#include "logger/logger.hpp"
#include "io/functors/average.hpp"
#include "plans-transforms/implemented_transformer.hpp"
#include "plans-solvers/implemented_equation.hpp"

#include "element.hpp"
#include "rezone.hpp"

namespace pisces
{
	enum initialize_element_flags {
		uniform_n = 0x01,
		uniform_m = 0x02
	};
	
	template <class datatype>
	class implemented_element : public pisces::element <datatype>
	{
	public:
		implemented_element (plans::axis i_axis_n, plans::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
		element <datatype> (i_name, 2, i_params, i_data, i_messenger_ptr, i_element_flags),
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
			init_timestep = i_params.get <datatype> ("time.init");
			mult_timestep = i_params.get <datatype> ("time.mult");
			timestep_safety = i_params.get <datatype> ("time.safety");
			next_timestep = 0.0;
			count = 0;
			previous = 0;
			
			element <datatype>::initialize (x_position, "x");
			element <datatype>::initialize (z_position, "z");
			if (i_params ["parallel.transform.subthreads"].IsDefined ()) {
				transform_threads = i_params.get <int> ("parallel.transform.subthreads");
			} else {
				transform_threads = 0;
			}
			
			grids [0] = data.get_grid (0);
			grids [1] = data.get_grid (1);
			
			for (typename data::data <datatype>::iterator iter = data.begin (); iter != data.end (); ++iter) {
				if ((iter->first != x_position) && (iter->first != z_position)) {
					DEBUG (ptr (iter->first));
					element <datatype>::add_solver (iter->first, std::shared_ptr <plans::equation <datatype> > (new plans::implemented_equation <datatype> (*grids [0], *grids [1], ptr (iter->first), &element_flags [state], &element_flags [iter->first])));
				
				}
			}
			
			TRACE ("Instantiated.");
		}
		
		virtual ~implemented_element () {}
		
		// /*!*******************************************************************
		//  * \brief Get the datatype reference to the named scalar
		//  *
		//  * \param name The integer name from the index enumeration
		//  *
		//  * \return A datatype reference to the first element of the named scalar
		//  *********************************************************************/
		// inline datatype& operator[] (int name) {
		// 	if (scalars.find (name) == scalars.end ()) {
		// 		FATAL ("Index " << name << " not found in element.");
		// 		throw 0;
		// 	}
		// 	return scalars [name] [0];
		// }
		
		inline datatype& operator() (int name, int i = 0, int j = 0) {
			return element <datatype>::operator() (name, i * m + j);
		}
		
		inline datatype* ptr (int name, int i = 0) {
			return element <datatype>::ptr (name, i);
		}

		inline datatype* ptr (int name, int i, int j) {
			return element <datatype>::ptr (name, i * m + j);
		}
	
		/*!*******************************************************************
		 * \copydoc element <datatype>::initialize ()
		 *********************************************************************/
		virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			// Size allowing for real FFT buffer
			// scalars [name].resize (grids [0]->get_ld () * m, 0.0);
			// if (name == x_position) {
			// 	for (int j = 0; j < m; ++j) {
			// 		linalg::copy (n, &((*grids [0]) [0]), ptr (name, 0, j), 1, m);
			// 	}
			// } else if (name == z_position) {
			// 	for (int i = 0; i < n; ++i) {
			// 		linalg::copy (m, &((*grids [1]) [0]), ptr (name, i, 0));
			// 	}
			// } else {
			// 	if (initial_conditions) {
			// 		if ((i_flags & uniform_m) && (i_flags & uniform_n)) {
			// 			for (int i = 0; i < n; ++i) {
			// 				for (int j = 0; j < m; ++j) {
			// 					(*this) (name, i, j) = *initial_conditions;
			// 				}
			// 			}
			// 		} else if (i_flags & uniform_m) {
			// 			for (int j = 0; j < m; ++j) {
			// 				linalg::copy (n, initial_conditions, ptr (name, 0, j), 1, m);
			// 			}
			// 		} else if (i_flags & uniform_n) {
			// 			for (int i = 0; i < n; ++i) {
			// 				linalg::copy (m, initial_conditions, ptr (name, i, 0));
			// 			}
			// 		} else {
			// 			linalg::copy (n * m, initial_conditions, ptr (name));
			// 		}
			// 	}
			// }
			//
			// if ((name != x_position) && (name != z_position)) {
			// 	element <datatype>::add_transform (name, std::shared_ptr <plans::transformer <datatype> > (new plans::implemented_transformer <datatype> (*grids [0], *grids [1], ptr (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal , &element_flags [state], &element_flags [name], transform_threads)));
			// }
			// if ((name != x_position) && (name != z_position)) {
			// 	DEBUG (ptr (name));
			// 	element <datatype>::add_solver (name, std::shared_ptr <plans::equation <datatype> > (new plans::implemented_equation <datatype> (*grids [0], *grids [1], ptr (name), &element_flags [state], &element_flags [name])));
			//
			// }
			TRACE ("Initialized.");
			return this->ptr (name);
		}
		
		virtual datatype calculate_timestep (int i, int j, io::formats::virtual_file *virtual_file = NULL) = 0;
		
		inline datatype calculate_min_timestep (io::formats::virtual_file *virtual_file = NULL) {
			double shared_min = max_timestep / timestep_safety;
			#pragma omp parallel 
			{
				double min = std::numeric_limits <double>::max ();
				#pragma omp for nowait
					for (int j = 1; j < (virtual_file ? virtual_file->dims ["z"] [1] : m) - 1; ++j) {
						for (int i = 0; i < (virtual_file ? virtual_file->dims ["z"] [0] : n); ++i) {
							min = std::min (calculate_timestep (i, j, virtual_file), min);
						}
					}
				#pragma omp critical 
				{
					shared_min = std::min (shared_min, min);
				}
			}
			shared_min *= timestep_safety;
			if (timestep == 0.0) {
				return init_timestep;
			}
			if (shared_min > mult_timestep * timestep) {
				if (next_timestep != 0.0 && std::min (next_timestep, shared_min) > mult_timestep * timestep) {
					next_timestep = 0.0;
					return std::min (mult_timestep * timestep, max_timestep);
				} else {
					next_timestep = shared_min;
					return timestep;
				}
			}
			next_timestep = 0.0;
			if (shared_min < timestep) {
				previous = count;
				return shared_min / mult_timestep;
			} else {
				return timestep;
			}
		}
		
		virtual io::formats::virtual_file *make_virtual_file (int flags = 0x00) {
			std::shared_ptr <io::formats::virtual_file> virtual_file (new io::formats::virtual_file);
			
			std::shared_ptr <io::output> virtual_output (new io::formatted_output <io::formats::two_d::virtual_format> (io::data_grid::two_d (n, m), "two_d/element/virtual_file", io::replace_file));
			data.setup_output (virtual_output, ::data::no_save);
			
			virtual_output->to_file ();
			return &io::virtual_files ["two_d/element/virtual_file"];
		}
		
		const int &get_mode () {
			return mode;
		}
		
		virtual std::shared_ptr <plans::grid <datatype>> generate_grid (plans::axis *axis, int index = -1) {
			if (index == 0) {
				return std::shared_ptr <plans::grid <datatype>> (new typename plans::horizontal::grid <datatype> (axis));
			} else if (index == 1 || index == -1) {
				return std::shared_ptr <plans::grid <datatype>> (new typename plans::vertical::grid <datatype> (axis));
			} else {
				throw 0;
			}
		}
		
		virtual io::formats::virtual_file *make_rezoned_virtual_file (datatype *positions, io::formats::virtual_file *old_virtual_file, int flags = 0x00) {
			plans::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <plans::grid <datatype>> vertical_grid = generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, old_virtual_file, &io::virtual_files ["two_d/element/new_virtual_file"]);
			
			return &io::virtual_files ["two_d/element/new_virtual_file"];
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
		using element <datatype>::scalars;
		using element <datatype>::grids;
		using element <datatype>::name;
		using element <datatype>::messenger_ptr;
		using element <datatype>::element_flags;
		using element <datatype>::timestep;
		using element <datatype>::axes;
		using element <datatype>::scalar_names;
		using element <datatype>::get_mode;
		using element <datatype>::duration;
		using element <datatype>::data;
		
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
		datatype init_timestep;
		datatype mult_timestep;
		datatype timestep_safety;
		datatype next_timestep;
		int count, previous;
		int transform_threads;
		static int mode;
	};
} /* pisces */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
