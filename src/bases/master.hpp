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
#include <algorithm>
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
		master (int i_id, int i_p, int i_n_elements, std::string parameter_filename) {
			// Read experiment parameters out of text file
			io::read_params_txt parameters (parameter_filename);
			inputParams = parameters.load_params();
			inputParams ["process_id"].asInt = i_id;
			inputParams ["total_processes"].asInt = i_p;

			// Inialize some experiment variables
			id = i_id;
			tsteps = inputParams ["timesteps"].asInt;
			
			/*
				TODO perhaps we can come up with a slightly more convenient method for retrieving from inputParams (e.g. inputParams.timesteps)?
			*/
			
			elements.resize (i_n_elements);
		}
		
		virtual ~master () {
			MTRACE ("Calling destructor.");
		}
		
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
			double t_timestep;
			std::vector <double> timesteps (inputParams ["total_processes"].asInt, 0.0);
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
				for (int j = 0; j < (int) elements.size (); ++j) {
					if (j == 0) {
						t_timestep = elements [j]->calculate_timestep ();
					} else {
						t_timestep = std::min (t_timestep, elements [j]->calculate_timestep ());
					}
				}
				if (inputParams ["total_processes"].asInt != 1) {
					MPI::COMM_WORLD.Gather (&t_timestep, 1, MPI_DOUBLE, &timesteps [0], 1, MPI_DOUBLE, 0);
					if (id == 0) {
						t_timestep = *std::min_element (timesteps.begin (), timesteps.end ());
					}
					MPI::COMM_WORLD.Bcast (&t_timestep, 1, MPI_DOUBLE, 0);
				}
				MTRACE ("Updating...");
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->update ();
				}
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->update_timestep (t_timestep);
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
