/*!***********************************************************************
 * \file master.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-05-14.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include <string>
#include "element.hpp"
#include "../utils/io.hpp"

namespace bases
{
	class master
	{
	public:
		master (int i_id, int i_n_elements, std::string parameter_filename) {
			// Read experiment parameters out of text file
			io::read_params_txt parameters (parameter_filename);
			inputParams = parameters.load_params();

			// Inialize some experiment variables
			id = i_id;
			tsteps = inputParams ["timesteps"].asInt;
			
			elements.resize (i_n_elements);
		}
		virtual ~master () {}
		
		virtual io::parameter_map& get_params () {
			return inputParams;
		}
		
		virtual void run () {
			MTRACE ("Beginning main...");
			for (int i = 0; i < tsteps; ++i) {
				MINFO ("Timestep " << i);
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->calculate ();
				}
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->send ();
				}
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->recv ();
				}
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->execute_boundaries ();
				}
				MTRACE ("Updating...");
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->update ();
				}
			}
		}

	protected:
		int id;
		int tsteps;
		
		std::vector <std::shared_ptr <element>> elements;
		io::parameter_map inputParams;
	};
} /* bases */
