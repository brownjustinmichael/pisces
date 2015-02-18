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

#include "exceptions.hpp"
#include "logger/logger.hpp"

namespace io
{
	/*!**********************************************************************
	 * \brief An extension of the YAML::Node class for convenience
	 * 
	 * This implementation uses the convenient shorthand of the '.' character delimiting nested YAML trees. This means that this character should be avoided in parameter names in your YAML file.
	 ************************************************************************/
	class parameters : public YAML::Node
	{
	public:
		// static parameters defaults;
		
		using YAML::Node::operator=; //!< Use the assignment operator from the super class
	
		/*!**********************************************************************
		 * \param file_name The parameter file from which the parameters should be loaded
		 ************************************************************************/
		parameters (std::string file_name = "", std::string defaults_file = DEFAULTS_FILE) {
			YAML::Node copy_node;
			YAML::Node::operator= (YAML::LoadFile (defaults_file));
			if (file_name != "") {
				copy_node = (YAML::LoadFile (file_name));
			}
			DEBUG ("HERE IS THE DEFAULT:" << (*this) ["output"]);
			YAML::Node::operator= (copy (copy_node, *this));
			DEBUG ("HERE IS THE FINAL:" << (*this) ["output"]);
		}
	
		virtual ~parameters () {}
		
		static YAML::Node copy (const YAML::Node &from, const YAML::Node to);
	
		/*!**********************************************************************
		 * \brief Overload the index operator for the YAML::Node
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * Note that this returns a copy of a YAML::Node object, which can be used to check and set parameters, but these objects do not treat parameter delimiters like this class.
		 * 
		 * \return A copy of the YAML::Node object at the given parameter.
		 ************************************************************************/
		YAML::Node operator[] (std::string key);
			
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
			if (operator[] (key).IsDefined ()) {
				return operator[] (key).as <datatype> ();
			} else {
				throw exceptions::key_does_not_exist (key);
			}
		}
	};
} /* io */

#endif /* end of include guard: PARAMETERS_HPP_A5EF2439 */
