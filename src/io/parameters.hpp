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
		using YAML::Node::operator=; //!< Use the assignment operator from the super class
	
		/*!**********************************************************************
		 * \param file_name The parameter file from which the parameters should be loaded
		 ************************************************************************/
		parameters (std::string file_name = "") {
			if (file_name != "") {
				YAML::Node::operator= (YAML::LoadFile (file_name));
			}
		}
	
		virtual ~parameters () {}
	
		/*!**********************************************************************
		 * \brief Overload the index operator for the YAML::Node
		 * 
		 * \param key The string representation of the parameter
		 * 
		 * Note that this returns a copy of a YAML::Node object, which can be used to check and set parameters, but these objects do not treat parameter delimiters like this class.
		 * 
		 * \return A copy of the YAML::Node object at the given parameter.
		 ************************************************************************/
		YAML::Node operator [] (std::string key) {
			std::istringstream ss (key);
			std::string token;
			std::getline (ss, token, '.');
			std::vector <YAML::Node> nodes;
			nodes.push_back (YAML::Node::operator [] (token));
			while (std::getline (ss, token, '.')) {
				nodes.push_back (nodes [(int) nodes.size () - 1] [token]);
			}
			return nodes [(int) nodes.size () - 1];
		}
			
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
		
		template <typename datatype>
		const datatype get (std::string key, datatype value) {
			if (operator[] (key).IsDefined ()) {
				return operator[] (key).as <datatype> ();
			} else {
				return value;
			}
		}
	};
} /* io */

#endif /* end of include guard: PARAMETERS_HPP_A5EF2439 */
