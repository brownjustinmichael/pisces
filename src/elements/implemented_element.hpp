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
	class implemented_element : public pisces::element
	{
	protected:
		using element::grids;
		using element::name;
		using element::messenger_ptr;
		using element::element_flags;
		using element::timestep;
		using element::axes;
		using element::get_mode;
		using element::duration;
		using element::data;
		
		int n; //!< The horizontal extent of the data
		int m; //!< The vertical extent of the data
		std::map <int, double> positions; //!< A vector of the edge positions
		std::map <int, int> excesses; //!< A vector of the edge positions
		std::vector<int> cell_n; //!< An integer array for tracking each cell number for output
		std::vector<int> cell_m; //!< An integer array for tracking each cell number for output
		std::vector <double> value_buffer_vec; //!< A vector for interpolating over the data in rezoning
		double* value_buffer; //!< A pointer to the value_buffer vector
		std::vector <double> inter_buffer_vec; //!< A vector for interpolating over the positions in rezoning
		double* inter_buffer; //!< A pointer to the inter_buffer vector
		
		double max_timestep; //!< The maximum timestep value
		double init_timestep; //!< The starting timestep value
		double mult_timestep; //!< The value by which to multiply or divide the timestep upon update
		double down_mult_timestep; //!< The value by which to multiply or divide the timestep upon update
		int check_every;
		double next_timestep; //!< The value of the next timestep
		int transform_threads; //!< The integer number of transform threads
		static int mode; //!< The mode of the simulation (from grids)
		
	public:
		typedef std::shared_ptr <element> (*element_function) (grids::axis, grids::axis, int, io::parameters&, data::data&, mpi::messenger*, int);

		static std::map <std::string, element_function> & registry ();

		/*!*******************************************************************
		* \param i_axis_n The axis for the horizontal direction
		* \param i_axis_m The axis for the vertical direction
		* \param i_name The string representation of the element
		* \param i_params The parameter object that contains the input parameters of the run
		* \param i_data An object that contains all the data in the simulation
		* \param i_messenger_ptr A pointer to a messenger object for inter-element communication
		* \param i_element_flags An integer set of global flags for the element
		*********************************************************************/
		implemented_element (grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) : 
		element (i_name, 2, i_params, i_data, i_messenger_ptr, i_element_flags),
		n (i_axis_n.get_n ()), m (i_axis_m.get_n ()) {
			TRACE ("Instantiating...");
			axes [0] = i_axis_n;
			axes [1] = i_axis_m;
			
			// Initialize the timestep
			max_timestep = i_params.get <double> ("time.max");
			init_timestep = i_params.get <double> ("time.init");
			mult_timestep = i_params.get <double> ("time.mult");
			down_mult_timestep = i_params.get <double> ("time.down_mult");
			increase_every = i_params.get <int> ("time.increase_every");
			next_timestep = 0.0;
			
			// Initialize x and z
			element::initialize ("x");
			element::initialize ("z");
			transform_threads = i_params.get <int> ("parallel.transform.subthreads");
			
			// Set up the grids
			grids [0] = std::shared_ptr <grids::grid> (new grids::horizontal::grid (&i_axis_n));
			grids [1] = std::shared_ptr <grids::grid> (new grids::vertical::grid (&i_axis_m));
			
			value_buffer_vec.resize (i_messenger_ptr->get_np () * m * n);
			value_buffer = &value_buffer_vec [0];
			inter_buffer_vec.resize (i_messenger_ptr->get_np () * m * n);
			inter_buffer = &inter_buffer_vec [0];
			
			// For every variable in the data object, add a corresponding equation and transformer
			for (data::data::iterator iter = data.begin (); iter != data.end (); ++iter) {
				if ((*iter != "x") && (*iter != "z")) {
					DEBUG ("Adding " << *iter);
					element::add_equation (*iter, std::shared_ptr <plans::solvers::equation > (new plans::solvers::implemented_equation (data [*iter], i_messenger_ptr)));
				}
			}
			
			TRACE ("Instantiated.");
		}
		
		virtual ~implemented_element () {}
		
		// /*!*******************************************************************
		//  * \brief Get the double reference to the named scalar
		//  *
		//  * \param name The integer name from the index enumeration
		//  *
		//  * \return A double reference to the first element of the named scalar
		//  *********************************************************************/
		// inline double& operator[] (int name) {
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
		inline double *operator() (std::string name, int i = 0, int j = 0) {
			return element::operator() (name, i * m + j);
		}
		
		/*!**********************************************************************
		 * \brief Get the pointer to the index of the specified variable
		 * 
		 * \param name The name of the variable to get
		 * \param i The index of the variable to get
		 * 
		 * \return The pointer to the index of the specified variable
		 ************************************************************************/
		inline double* ptr (std::string name, int i = 0) {
			return element::ptr (name, i);
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
		inline double* ptr (std::string name, int i, int j) {
			return element::ptr (name, i * m + j);
		}
		
		/*!**********************************************************************
		 * @brief Calculate the timestep for a given cell
		 * 
		 * @param i The horizontal index
		 * @param j The vertical index
		 * @param virtual_file A pointer to a virtual file for which to calculate the timestep, if NULL use the current state
		 * @return The timestep for the cell
		 ************************************************************************/
		virtual double calculate_timestep (int i, int j, formats::virtual_file *virtual_file = NULL) = 0;
		
		/*!**********************************************************************
		 * \copydoc element::calculate_min_timestep
		 ************************************************************************/
		inline double calculate_min_timestep (formats::virtual_file *virtual_file = NULL, bool limiters = true) {
			double shared_min = max_timestep;
			static int count = 0;
			
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
					return std::min (init_timestep, max_timestep);
				}
				if (shared_min * down_mult_timestep > timestep) {
					// If the minimum is larger than the current, increase the timestep
					count++;
					if (count % increase_every == 0) {
						return std::min (mult_timestep * timestep, max_timestep);
					} else {
						return timestep;
					}
				}
				if (shared_min < timestep) {
					// If the minimum is lower than the current, decrease the timestep
					count = 0;
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
		virtual std::shared_ptr <grids::grid> generate_grid (grids::axis *axis_ptr, int index = -1) {
			if (index == 0) {
				return std::shared_ptr <grids::grid> (new grids::horizontal::grid (axis_ptr));
			} else if (index == 1 || index == -1) {
				return std::shared_ptr <grids::grid> (new grids::vertical::grid (axis_ptr));
			} else {
				throw 0;
			}
		}
		
		/*!**********************************************************************
		 * \copydoc element::make_rezoned_virtual_file
		 ************************************************************************/
		virtual formats::virtual_file *make_rezoned_virtual_file (double *positions, formats::virtual_file *virtual_file_ptr, int flags = 0x00) {
			grids::axis vertical_axis (m, positions [messenger_ptr->get_id ()], positions [messenger_ptr->get_id () + 1], messenger_ptr->get_id () == 0 ? 0 : 1, messenger_ptr->get_id () == messenger_ptr->get_np () - 1 ? 0 : 1);
			std::shared_ptr <grids::grid> vertical_grid = generate_grid (&vertical_axis);
			
			pisces::rezone (messenger_ptr, &*(grids [1]), &*vertical_grid, virtual_file_ptr, &formats::virtual_files ["two_d/element/new_virtual_file"], value_buffer, inter_buffer);
			
			return &formats::virtual_files ["two_d/element/new_virtual_file"];
		}
		
		/*
			TODO Combine these?
		*/
		
		/*!**********************************************************************
		 * \copydoc element::get_zoning_positions
		 ************************************************************************/
		virtual void get_zoning_positions (double *positions) {
			if (messenger_ptr->get_id () == 0) {
				double temp [messenger_ptr->get_np () * 2];
				temp [0] = axes [1].get_position_0 ();
				temp [1] = axes [1].get_position_n ();
				messenger_ptr->gather <double> (2, temp, temp);
				for (int i = 0; i < messenger_ptr->get_np (); ++i) {
					positions [i] = temp [2 * i];
					DEBUG ("REZONING " << positions [i]);
				}
				positions [messenger_ptr->get_np ()] = temp [messenger_ptr->get_np () * 2 - 1];
				DEBUG ("REZONING " << positions [messenger_ptr->get_np ()]);
				messenger_ptr->bcast <double> (messenger_ptr->get_np () + 1, positions);
			} else {
				double temp [2] = {axes [1].get_position_0 (), axes [1].get_position_n ()};
				messenger_ptr->gather <double> (2, temp, temp);
				messenger_ptr->bcast <double> (messenger_ptr->get_np () + 1, positions);
			}
		}

		static void registrate (std::string const & name, element_function fp) {
			registry () [name] = fp;
		}

		static std::shared_ptr <element> instance (std::string const & name, grids::axis i_axis_n, grids::axis i_axis_m, int i_name, io::parameters& i_params, data::data &i_data, mpi::messenger* i_messenger_ptr, int i_element_flags) {
			auto it = registry ().find (name);
			element_function fp = (it == registry ().end () ? NULL : (it->second));
			if (fp) {
				return fp (i_axis_n, i_axis_m, i_name, i_params, i_data, i_messenger_ptr, i_element_flags);
			} else {
				return std::shared_ptr <element> ();
			}
		}

		template <typename type>
		struct registrar
		{
			explicit registrar (std::string const & name) {
				implemented_element::registrate (name, &type::instance);
			}
		};
	};
} /* pisces */

#endif /* end of include guard: ELEMENT_TWO_D_HPP_CJ68F4IB */
