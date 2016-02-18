/*!**********************************************************************
 * \file implemented_data.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-09.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef IMPLEMENTED_DATA_HPP_574F2B9E
#define IMPLEMENTED_DATA_HPP_574F2B9E

#include <string>
#include <map>
#include <memory>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <unistd.h>

#include "linalg/utils.hpp"
#include "data.hpp"

#include "io/functors/slice.hpp"
#include "io/functors/product.hpp"
#include "io/input.hpp"
#include "io/formats/format.hpp"
#include "io/functors/functor.hpp"
#include "io/output.hpp"
#include "io/functors/average.hpp"
#include "io/formats/exceptions.hpp"
#include "io/formats/netcdf.hpp"
#include "plans/grids/grid.hpp"
#include "plans-transforms/implemented_transformer.hpp"

namespace data
{
	/*!**********************************************************************
	 * \brief An implemeted form of data in 2D
	 ************************************************************************/
	class implemented_data : public data
	{
	protected:
		using data::iterator;
		using data::weights;
		using data::variables;
		using data::flags;
		
		std::shared_ptr <grids::grid> grid_n; //!< The horizontal grid object
		std::shared_ptr <grids::grid> grid_m; //!< The vertical grid object
		int n; //!< The horizontal extent of the data
		int m; //!< The vertical extent of the data

		std::vector <double> area; //!< A calculation of the area encompassed by each cell in the grid
		
	public:
		using data::transformers;
		using data::duration;

		typedef std::shared_ptr <data> (*data_function) (grids::axis*, grids::axis*, int, int, io::parameters&);

		static std::map <std::string, data_function> & registry ();
		
		/*!**********************************************************************
		 * \param i_axis_n The horizontal axis object
		 * \param i_axis_m The vertical axis object
		 * @param i_params The io::parameters object associated with the run
		 * \param i_name The integer name of the element
		 * \param dump_file The string name of the dump file (no extension)
		 * \param dump_directory The string directory to store the dump
		 * \param dump_every The frequency to dump to file
		 ************************************************************************/
		implemented_data (grids::axis *i_axis_n, grids::axis *i_axis_m, io::parameters &i_params, int i_name = 0, std::string dump_file = "", std::string dump_directory = "./", int dump_every = 1);
		
		virtual ~implemented_data () {}

		static void registrate (std::string const & name, data_function fp) {
			registry () [name] = fp;
		}

		static std::shared_ptr <data> instance (std::string const & name, grids::axis* i_axis_n, grids::axis* i_axis_m, int i_name, int n_elements, io::parameters& i_params);

		template <typename type>
		struct registrar
		{
			explicit registrar (std::string const & name) {
				implemented_data::registrate (name, &type::instance);
			}
		};

		/*!**********************************************************************
		 * \copydoc data::setup_profile
		 ************************************************************************/
		virtual void setup_profile (std::shared_ptr <io::output> output_ptr, int flags = 0x00);

		/**
		 * @brief Set up the data from a YAML Node
		 * @details This calls the setup_from method of the superclass, but here the data_grid object is constructed for you for convenience
		 * 
		 * @param input_params The YAML Node to construct the input stream from
		 * @return A shared pointer to the resulting input stream object
		 */
		template <class format>
		std::shared_ptr <io::input> setup_from (YAML::Node input_params) {
			const formats::data_grid grid = formats::data_grid::two_d (n, m);

			return data::template setup_from <format> (input_params, grid);
		}

		/**
		 * @brief Set up an output stream from a YAML Node
		 * @details This calls the setup_output_from method of the superclass, but here the data_grid object is constructed for you for convenience
		 * 
		 * @param output_params The YAML Node from which the output stream should be constructed
		 * @param state The state from which anything in the output stream should be read
		 * @param flags Any flags to the output stream setup that are desired to be passed
		 * @return A shared pointer to the newly constructed output stream 
		 */
		template <class format>
		std::shared_ptr <io::output> setup_output_from (YAML::Node output_params, int state = real_real, int flags = 0x00) {
			const formats::data_grid grid = formats::data_grid::two_d (n, m);

			return data::template setup_output_from <format> (output_params, grid, state, flags);
		}

		/**
		 * @brief Construct a functor that returns the max of a given string variable
		 * 
		 * @param variable The string variable from which the max should be gathered
		 * @return The max functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_max (std::string variable);

		/**
		 * @brief Construct a functor that returns the average of a given string variable
		 * 
		 * @param variable The string variable from which the average should be calculated
		 * @return The average functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_avg (std::string variable);

		/**
		 * @brief Construct a functor that takes the average of the derivative for a given string variable in the vertical direction
		 * 
		 * @param variable The string variable from which the average of the derivative should be calculated
		 * @return The derivative functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_deriv (std::string variable, int slice_index = -1);

		/**
		 * @brief Construct a flux functor for a given string variable
		 * 
		 * @param variable The string variable for which to construct the flux
		 * @param velocity The string representation of the velocity variable needed to construct the flux
		 * 
		 * @return The flux functor, which should be added to an output stream.
		 */
		std::shared_ptr <functors::functor> output_flux (std::string variable, std::string velocity, int slice_index = -1);

		/**
		 * @brief Construct a flux functor for a given string variable
		 * 
		 * @param variable The string variable for which to construct the flux
		 * @param velocity The string representation of the velocity variable needed to construct the flux
		 * 
		 * @return The flux functor, which should be added to an output stream.
		 */
		std::shared_ptr <functors::functor> output_avg_flux (std::string variable, std::string velocity);

		/**
		 * @brief Construct a functor that returns the profile of a given string variable
		 * 
		 * @param variable The string variable from which the profile should be calculated
		 * @return The profile functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_prof (std::string variable);

		/**
		 * @brief Construct a functor that returns the profile of the vertical flux of a given string variable
		 * 
		 * @param variable The string variable from which the profile should be calculated
		 * @param velocity The string representation of the velocity variable needed to construct the flux
		 * @return The profile functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_flux_prof (std::string variable, std::string velocity);

		/**
		 * @brief Construct a functor that returns the profile of the derivative of a given string variable
		 * 
		 * @param variable The string variable from which the profile should be calculated
		 * @return The profile functor, which should be added to an output stream
		 */
		std::shared_ptr <functors::functor> output_deriv_prof (std::string variable);

		/**
		 * @brief Index the data in two dimensions and return the appropriate pointer
		 * 
		 * @param name The name of the variable to index
		 * @param state The state of the variable to index
		 * @param i The x index of the variable to index
		 * @param j The z index of the variable to index
		 * @return A pointer to the data cell for the given variable and state
		 */
		double *operator() (const std::string &name, int state = 0, int i = 0, int j = 0);

	protected:
		/*!**********************************************************************
		 * \copydoc data::_initialize
		 ************************************************************************/
		virtual grids::variable &_initialize (std::string name, int i_flags = 0x00);
	};	
}

#endif /* end of include guard: IMPLEMENTED_DATA_HPP_574F2B9E */
