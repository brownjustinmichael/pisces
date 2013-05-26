/*!***********************************************************************
 * \file master.hpp
 * Spectral Element
 * 
 * This file contains the master class, which contains a set of elements 
 * and slave processes. A single run of the code may in fact have multiple 
 * masters, each with a number of elements and/or slaves.
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
	/*!*******************************************************************
	 * \brief The class that sets up and runs elements
	 * 
	 * The master class contains a set of elements and slave processes. 
	 * Generally, the master will take care of the administrative details
	 * where the slaves will mostly just run linear algebra.
	 *********************************************************************/
	class master
	{
	public:
		/*!*******************************************************************
		 * \param i_id The integer process id
		 * \param i_n_elements The number of elements in the process
		 * \param parameter_filename The string location of the parameter file
		 *********************************************************************/
		master (int i_id, int i_n_elements, std::string parameter_filename) {
			// Read experiment parameters out of text file
			io::read_params_txt parameters (parameter_filename);
			inputParams = parameters.load_params();

			// Inialize some experiment variables
			id = i_id;
			tsteps = inputParams ["timesteps"].asInt;
			
			/*
				TODO perhaps we can come up with a slightly more convenient method for retrieving from inputParams (e.g. inputParams.timesteps)?
			*/
			
			elements.resize (i_n_elements);
		}
		
		virtual ~master () {}
		
		/*!*******************************************************************
		 * \brief Get the parameter map object associated with the master object
		 * 
		 * \returns The parameter map object
		 *********************************************************************/
		virtual io::parameter_map& get_params () {
			return inputParams;
		}
		
		/*!*******************************************************************
		 * \brief Run the simulation
		 * 
		 * Go through each element in order and execute each.
		 *********************************************************************/
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
		int id; //!< The integer processor id
		int tsteps; //!< The integer total number of timesteps to take
		
		std::vector <std::shared_ptr <element>> elements; //!< A vector containing shared pointers to the contained elements
		io::parameter_map inputParams; //!< The parameter map object containing the input parameters
	};
} /* bases */
