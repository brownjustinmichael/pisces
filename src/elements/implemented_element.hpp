/*!**********************************************************************
 * \file implemented_element.hpp
 * /Users/justinbrown/Dropbox/pisces/src
 * 
 * Created by Justin Brown on 2013-08-21.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_TWO_D_HPP_CJ68F4IB
#define ELEMENT_TWO_D_HPP_CJ68F4IB

#include "mpi/messenger.hpp"

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
	/*!**********************************************************************
	 * \brief A set of flags to be used when initializing element data
	 ************************************************************************/
	enum initialize_element_flags {
		uniform_n = 0x01,
		uniform_m = 0x02
	};
	
	/*!**********************************************************************
	 * \brief An implementation of the element class in 2D
	 ************************************************************************/
	template <class datatype>
	class implemented_element : public pisces::element <datatype>
	{
	protected:
		using element <datatype>::grids;
		using element <datatype>::name;
		using element <datatype>::messenger_ptr;
		using element <datatype>::element_flags;
		using element <datatype>::timestep;
		using element <datatype>::axes;
		using element <datatype>::get_mode;
		using element <datatype>::duration;
		using element <datatype>::data;
		
		int n; //!< The horizontal extent of the data
		int m; //!< The vertical extent of the data
		std::map <int, datatype> positions; //!< A vector of the edge positions
		std::map <int, int> excesses; //!< A vector of the edge positions
		std::vector<int> cell_n; //!< An integer array for tracking each cell number for output
		std::vector<int> cell_m; //!< An integer array for tracking each cell number for output
		std::vector <datatype> value_buffer_vec;
		datatype* value_buffer;
		std::vector <datatype> inter_buffer_vec;
		datatype* inter_buffer;
		
		datatype max_timestep; //!< The maximum timestep value
		datatype init_timestep; //!< The starting timestep value
		datatype mult_timestep; //!< The value by which to multiply or divide the timestep upon update
		datatype down_mult_timestep; //!< The value by which to multiply or divide the timestep upon update
		datatype next_timestep; //!< The value of the next timestep
		int transform_threads; //!< The integer number of transform threads
		static int mode; //!< The mode of the simulation (from grids)
		
	public:
		/*!*******************************************************************
		* \param i_axis_n The axis for the horizontal direction
		* \param i_axis_m The axis for the vertical direction
		* \param i_name The string representation of the element
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_data An object that contains all the data in the simulation
		* \param i_messenger_ptr A pointer to a messenger object for inter-element communication
		* \param i_element_flags An integer set of global flags for the element
		*********************************************************************/
		implemented_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data <datatype> &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
		element <datatype> (i_name, 2, i_params, i_data, i_messenger_ptr, i_element_flags),
		n (i_axis_n.get_n ()), m (i_axis_m.get_n ()) {
			TRACE ("Instantiating...");
			axes [0] = i_axis_n;
			axes [1] = i_axis_m;
			
			// Initialize the timestep
			max_timestep = i_params.get <datatype> ("time.max");
			init_timestep = i_params.get <datatype> ("time.init");
			mult_timestep = i_params.get <datatype> ("time.mult");
			down_mult_timestep = i_params.get <datatype> ("time.down_mult");
			next_timestep = 0.0;
			
			// Initialize x and z
			element <datatype>::initialize ("x");
			element <datatype>::initialize ("z");
			transform_threads = i_params.get <int> ("parallel.transform.subthreads");
			
			// Set up the grids
			grids [0] = std::shared_ptr <grids::grid <datatype>> (new typename grids::horizontal::grid <datatype> (&i_axis_n));
			grids [1] = std::shared_ptr <grids::grid <datatype>> (new typename grids::vertical::grid <datatype> (&i_axis_m));
			
			value_buffer_vec.resize (i_messenger_ptr->get_np () * m * n);
			value_buffer = &value_buffer_vec [0];
			inter_buffer_vec.resize (i_messenger_ptr->get_np () * m * n);
			inter_buffer = &inter_buffer_vec [0];
			
			// For every variable in the data object, add a corresponding equation and transformer
			for (typename data::data <datatype>::iterator iter = data.begin (); iter != data.end (); ++iter) {
				if ((*iter != "x") && (*iter != "z")) {
					DEBUG ("Adding " << *iter);
					element <datatype>::add_equation (*iter, std::shared_ptr <plans::solvers::equation <datatype> > (new plans::solvers::implemented_equation <datatype> (data [*iter], &element_flags ["element"], &element_flags [*iter], i_messenger_ptr)));
					element <datatype>::transforms.push_back (*iter);
					element <datatype>::transformers [*iter] = std::shared_ptr <plans::transforms::transformer <datatype> > (new plans::transforms::implemented_transformer <datatype> (*grids [0], *grids [1], data (*iter), NULL, plans::transforms::forward_vertical | plans::transforms::forward_horizontal | plans::transforms::inverse_vertical | plans::transforms::inverse_horizontal , &(data.flags ["element"]), &(data.flags [*iter]), transform_threads));
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
		
		/*!**********************************************************************
		 * \copydoc implemented_element::ptr
		 * 
		 * \param j The vertical index of the variable to get
		 ************************************************************************/
		inline datatype *operator() (std::string name, int i = 0, int j = 0) {
			return element <datatype>::operator() (name, i * m + j);
		}
		
		/*!**********************************************************************
		 * \brief Get the pointer to the index of the specified variable
		 * 
		 * \param name The name of the variable to get
		 * \param i The index of the variable to get
		 * 
		 * \return The pointer to the index of the specified variable
		 ************************************************************************/
		inline datatype* ptr (std::string name, int i = 0) {
			return element <datatype>::ptr (name, i);
		}
		
		/*!**********************************************************************
		 * \brief Get the pointer to the index of the specified variable
		 * 
		 * \param name The name of the variable to get
		 * \param i The horizontal index of the variable to get
		 * \param j The vertical index of the variable to get
		 * 
		 * \return The pointer to the index of the specified variable
		 ************************************************************************/
		inline datatype* ptr (std::string name, int i, int j) {
			return element <datatype>::ptr (name, i * m + j);
		}
		
		/*!**********************************************************************
		 * \copydoc element::calculate_timestep
		 ************************************************************************/
		virtual datatype calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL) = 0;
		
		/*!**********************************************************************
		 * \copydoc element::calculate_min_timestep
		 ************************************************************************/
		inline datatype calculate_min_timestep (formats::virtual_file *virtual_file = NULL, bool limiters = true) {
			double shared_min = max_timestep;
			
			// Calculate the minimum timestep in parallel
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
			
			// Using the minimum timestep, check whether the code should increase or reduce the timestep
			if (limiters) {
				DEBUG ("Desired: " << shared_min << " Current: " << timestep);
				if (timestep == 0.0) {
					return init_timestep;
				}
				if (shared_min * down_mult_timestep > timestep) {
					// If the minimum is larger than the current, increase the timestep
					next_timestep = 0.0;
					return std::min (mult_timestep * timestep, max_timestep);
				}
				if (shared_min < timestep) {
					// If the minimum is lower than the current, decrease the timestep
					return shared_min * down_mult_timestep;
				} else {
					return timestep;
				}
			} else {
				return shared_min;
			}
		}
		
		/*!**********************************************************************
		 * \copydoc element::make_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_virtual_file (int flags = 0x00) {
			std::shared_ptr <formats::virtual_file> virtual_file (new formats::virtual_file);
			
			std::shared_ptr <io::output> virtual_output (new io::formatted_output <formats::virtual_format> (formats::data_grid::two_d (n, m), "two_d/element/virtual_file", formats::replace_file));
			data.setup_output (virtual_output, ::data::no_save);
			
			virtual_output->to_file ();
			return &formats::virtual_files ["two_d/element/virtual_file"];
		}
		
		/*!**********************************************************************
		 * \copydoc element::get_mode
		 ************************************************************************/
		const int &get_mode () {
			return mode;
		}
		
		/*!**********************************************************************
		 * \copydoc element::generate_grid
		 ************************************************************************/
		virtual std::shared_ptr <grids::grid <datatype>> generate_grid (grids::axis *axis_ptr, int index = -1) {
			if (index == 0) {
				return std::shared_ptr <grids::grid <datatype>> (new typename grids::horizontal::grid <datatype> (axis_ptr));
			} else if (index == 1 || index == -1) {
				return std::shared_ptr <grids::grid <datatype>> (new typename grids::vertical::grid <datatype> (axis_ptr));
			} else {
				throw 0;
			}
		}
		
		/*!**********************************************************************
		 * \copydoc element::make_rezoned_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_rezoned_virtual_file (datatype *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00) {
			grids::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <grids::grid <datatype>> vertical_grid = generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, virtual_file_ptr, &formats::virtual_files ["two_d/element/new_virtual_file"], value_buffer, inter_buffer);
			
			return &formats::virtual_files ["two_d/element/new_virtual_file"];
		}
		
		/*
			TODO Combine these?
		*/
		
		/*!**********************************************************************
		 * \copydoc element::get_zoning_positions
		 ************************************************************************/
		virtual void get_zoning_positions (datatype *positions) {
			if (messenger_ptr->get_id () == 0) {
				datatype temp [messenger_ptr->get_np () * 2];
				temp [0] = axes [1].get_position_0 ();
				temp [1] = axes [1].get_position_n ();
				messenger_ptr->template gather <datatype> (2, temp, temp);
				for (int i = 0; i < messenger_ptr->get_np (); ++i) {
					positions [i] = temp [2 * i];
					DEBUG ("REZONING " << positions [i]);
				}
				positions [messenger_ptr->get_np ()] = temp [messenger_ptr->get_np () * 2 - 1];
				DEBUG ("REZONING " << positions [messenger_ptr->get_np ()]);
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			} else {
				datatype temp [2] = {axes [1].get_position_0 (), axes [1].get_position_n ()};
				messenger_ptr->template gather <datatype> (2, temp, temp);
				messenger_ptr->template bcast <datatype> (messenger_ptr->get_np () + 1, positions);
			}
		}
	};
} /* pisces */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
