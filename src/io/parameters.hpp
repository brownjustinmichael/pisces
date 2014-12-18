/*!**********************************************************************
 * \file parameters.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PARAMETERS_HPP_A5EF2439
#define PARAMETERS_HPP_A5EF2439

#include <fstream>

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
		/*!**********************************************************************
		 * \param file_name The parameter file from which the parameters should be loaded
		 ************************************************************************/
		parameters (std::string file_name = "") {
			if (file_name != "") {
			    std::ifstream fin (file_name);
			    YAML::Parser parser(fin);
				
			    parser.GetNextDocument (*this);
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
		const YAML::Node &operator [] (std::string key) {
			std::istringstream ss (key);
			std::string token;
			std::getline (ss, token, '.');
			const YAML::Node *node = &(YAML::Node::operator [] (token));
			while (std::getline (ss, token, '.')) {
				node = &((*node) [token]);
			}
			return *node;
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
			datatype value;
			(*this) [key] >> value;
			return value;
		}
		
		template <typename datatype>
		const datatype get (std::string key, datatype default_value) {
			datatype value;
			try {
				(*this) [key] >> value;
			} catch (YAML::RepresentationException &e) {
				return value;
			}
			return default_value;
		}
	};
} /* io */

#endif /* end of include guard: PARAMETERS_HPP_A5EF2439 */
