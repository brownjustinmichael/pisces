/*!**********************************************************************
 * \file parameters.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PARAMETERS_HPP_A5EF2439
#define PARAMETERS_HPP_A5EF2439

#include <yaml-cpp/yaml.h>

#include "versions/version.hpp"
#include "exceptions.hpp"
#include "logger/logger.hpp"

#ifndef DEFAULTS_FILE
#define DEFAULTS_FILE ""
#endif

namespace io
{
	/*!**********************************************************************
	 * \brief An extension of the YAML::Node class for convenience
	 * 
	 * This implementation uses the convenient shorthand of the '.' character delimiting nested YAML trees. This means that this character should be avoided in parameter names in your YAML file. Note that these can be changed at any time.
	 ************************************************************************/
	class parameters : public YAML::Node
	{
	public:
		using YAML::Node::operator=; // Use the assignment operator from the super class
	
		/*!**********************************************************************
		 * \brief An alias class designed to be a reference to a parameter node
		 * 
		 * The alias class is useful when iterating through keys in the parameter class. Since such iterations will store YAML::Node objects, you lose the ability to use '.' separated keys and any knowledge of the parameters on the whole.
		 ************************************************************************/
		class alias
		{
		private:
			parameters &params; //!< The parameters object to alias
			std::string key; //!< The key of params at which the alias is placed
			
		public:
			/*!**********************************************************************
			 * \param i_params The parameters object to alias
			 * \param i_key The key of the parameters object at which to make the alias
			 ************************************************************************/
			alias (parameters &i_params, std::string i_key) : params (i_params), key (i_key) {}
			
			virtual ~alias () {}
			
			YAML::Node operator [] (std::string i_key) {
				return params [key + "." + i_key];
			}
		};
	
		/*!**********************************************************************
		 * \param file_name The parameter file from which the parameters should be loaded
		 ************************************************************************/
		parameters (std::string file_name = "", std::string defaults_file = DEFAULTS_FILE) {
			TRACE ("Constructing parameters");
			YAML::Node copy_node;
			
			// Load the defaults file (by default, the DEFAULTS_FILE macro)
			if (defaults_file != "") {
				YAML::Node::operator= (YAML::LoadFile (defaults_file));
			}
			
			// If a file_name is specified, load it into copy
			if (file_name != "") {
				copy_node = (YAML::LoadFile (file_name));
			}
			
			// Override anything specified in file_name
			YAML::Node::operator= (copy (copy_node, *this));
			TRACE ("Constructed");
		}
	
		virtual ~parameters () {}
		
		/*!**********************************************************************
		 * \brief Get the version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.1.0");
			return version;
		}
		
		/*!**********************************************************************
		 * \brief A static method that overrides every value in to with those from from and returns a new node
		 * 
		 * \param from The YAML::Node object from which to copy the information
		 * \param to The YAML::Node object to which to copy the information
		 * 
		 * For each variable present in from, that variable will be created or overwritten in to. This is done recursively in order to copy maps from node to node.
		 * 
		 * \return A new node that is a copy of to with overwritted components of from
		 ************************************************************************/
		static YAML::Node copy (const YAML::Node &from, const YAML::Node to);
	
		/*!**********************************************************************
		 * \brief Overload the index operator for the YAML::Node
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * Note that this returns a copy of a YAML::Node object, which can be used to check and set parameters, but these objects do not treat parameter delimiters like this class. This operation will also search for linked keys. For example, "equations.x_velocity" has no "diffusion" key, but it does have a "link" key that points to "equations.velocity". That key does have a "diffusion" key, so "equations.x_velocity.diffusion" will return the value of "equations.velocity.diffusion".
		 * 
		 * \return A copy of the YAML::Node object at the given parameter.
		 ************************************************************************/
		YAML::Node operator[] (std::string key) const;
			
		/*!**********************************************************************
		 * \brief Get the parameter associated with the given key
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * This method is included for convenience. It avoids some of the messy YAML syntax and takes advantage of the overwritted index operator. This method throws a key_does_not_exist exception if it cannot find the key.
		 * 
		 * \return Returns the value of the given parameter
		 ************************************************************************/
		template <typename datatype>
		const datatype get (std::string key) {
			YAML::Node node = operator[] (key);
			if (operator[] (key).IsDefined ()) {
				return operator[] (key).as <datatype> ();
			} else {
				ERROR ("Key " << key << " does not exist.");
				throw exceptions::key_does_not_exist (key);
			}
		}
	};
} /* io */

#endif /* end of include guard: PARAMETERS_HPP_A5EF2439 */
