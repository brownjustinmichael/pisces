/*!**********************************************************************
 * \file parameters.cpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2015-01-27.
 * Copyright 2015 Justin Brown. All rights reserved.
 ************************************************************************/

#include "parameters.hpp"
#include <iostream>

namespace io
{
	YAML::Node parameters::copy (const YAML::Node &from, const YAML::Node to) {
		/*
			TODO This currently doesn't allow for sequences
		*/
		YAML::Node copyNode;
		
		// If the value to copy is a scalar, return it; it should overwrite whatever is present in to
		if (from.IsScalar ()) {
			return from;
		}
		
		/*
			TODO This is not the desired behavior: a map in from will not overwrite a scalar in to...
		*/
		if (to.IsScalar ()) {
			return to;
		}
		
		// Make a copy of to containing the entire tree
		if (to.IsMap ()) {
			for (YAML::const_iterator iter = to.begin (); iter != to.end (); ++iter) {
				copyNode [(iter->first).as <std::string> ()] = iter->second;
			}
		}
		
		// For every key value pair in from
		if (from.IsMap ()) {
			for (YAML::const_iterator iter = from.begin (); iter != from.end (); ++iter) {
				// If the corresponding key exists in to, call copy to correctly check for overwriting sub maps; otherwise, copy it
				if (copyNode [(iter->first).as <std::string> ()].IsDefined ()) {
					copyNode [(iter->first).as <std::string> ()] = copy (iter->second, copyNode [(iter->first).as <std::string> ()]);
				} else {
					copyNode [(iter->first).as <std::string> ()] = iter->second;
				}
			}
		}
		return copyNode;
	}
	
	YAML::Node parameters::operator[] (std::string key) const {
		DEBUG ("Looking for " << key);

		// Tokenize the incoming string key using '.' as a delimeter
		std::istringstream ss (key);
		std::string token;
		std::getline (ss, token, '.');
		std::vector <std::string> tokens;
		std::vector <YAML::Node> nodes;
		
		// Iterate through the tokens, grabbing the relevant nodes
		tokens.push_back (token);
		nodes.push_back (YAML::Node::operator[] (token));
		while (std::getline (ss, token, '.')) {
			tokens.push_back (token);
			nodes.push_back (nodes [(int) nodes.size () - 1] [token]);
		}

		if (tokens.size () == 1) {
			return nodes [0];
		}
		
		// Check that the final node is actually defined
		int i = (int) nodes.size () - 1;
		std::string inner_tokens = "";
		if (!nodes [i].IsDefined ()) {
			// If not, search back through the nodes for a "link" key
			for (int j = i - 1; j >= 0; --j) {
				inner_tokens = tokens [j + 1] + inner_tokens;
				if (nodes [j].IsDefined () && nodes [j] ["link"].IsDefined ()) {
					// If a "link" key is found, check whether that key has the appropriate key
					if (parameters::operator[] (nodes [j] ["link"].as <std::string> () + "." + inner_tokens).IsDefined ()) {
						return parameters::operator [] (nodes [j] ["link"].as <std::string> () + "." + inner_tokens);
					}
				}
				inner_tokens = "." + inner_tokens;
			}
		}
		return nodes [i];
	}
} /* io */
