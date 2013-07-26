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
#include <iomanip>
#include <algorithm>
#include "element.hpp"
#include "../utils/io.hpp"
#include "../utils/messenger.hpp"

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
		master (int i_id, int i_p, int i_n_elements, std::string parameter_filename, utils::messenger* i_messenger_ptr) {
			// Read experiment parameters out of text file
			io::read_params_txt parameters (parameter_filename);
			messenger_ptr = i_messenger_ptr;
			inputParams = parameters.load_params();
			inputParams ["process_id"].asInt = i_id;
			inputParams ["total_processes"].asInt = i_p;

			// Inialize some experiment variables
			id = i_id;
			tsteps = inputParams ["timesteps"].asInt;
			count = 0;
			
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
			int info;
			int status;
			std::vector <int> stati ((int) elements.size (), 0);
			std::vector <double> timesteps (inputParams ["total_processes"].asInt, 0.0);
			for (int i = 0; i < tsteps; ++i) {
				MINFO ("Timestep " << i);
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->calculate ();
					elements [j]->output ();
					elements [j]->execute_boundaries ();
					if (j == 0) {
						t_timestep = elements [j]->calculate_timestep ();
					} else {
						t_timestep = std::min (t_timestep, elements [j]->calculate_timestep ());
					}
				}
				messenger_ptr->min (&t_timestep);
				MTRACE ("Updating...");
				status = 0;
				for (int j = 0; j < (int) elements.size (); ++j) {
					stati [j] = 0;
				}
				utils::scale (count * count, 0.0, &global_matrix [0]);
				for (int j = 0; j < count; ++j) {
					ipiv [j] = 0;
				}
				while (status != 1) {
					for (int j = 0; j < (int) elements.size (); ++j) {
						elements [j]->update_globals (count, &global_matrix [0], &global_rhs [0], &stati [j]);
						if (status == 0) {
							if (j == 0 && stati [0] == 1) {
								status = 1;
							}
						} else if (status == 1 && stati [j] != 1) {
							status = 0;
						}
					}
				}
				utils::matrix_factorize (count, count, &global_matrix [0], &ipiv [0], &info);
				utils::matrix_solve (count, &global_matrix [0], &ipiv [0], &global_rhs [0], &info);
				for (int j = 0; j < (int) elements.size (); ++j) {
					elements [j]->update_from_globals (&global_rhs [0]);
				}
				// for (int k = 0; k < 2; ++k) {
				// 	for (int j = 0; j < (int) elements.size (); ++j) {
				// 		elements [j]->attempt_update ();
				// 		elements [j]->calculate_bounds ();
				// 		elements [j]->send_bounds ();
				// 	}
				// 	for (int j = 0; j < (int) elements.size (); ++j) {
				// 		elements [j]->recv_bounds ();
				// 		elements [j]->calculate_error ();
				// 	}
				// 	for (int j = 0; j < (int) elements.size (); ++j) {
				// 		elements [j]->send_error ();
				// 	}
				// 	for (int j = 0; j < (int) elements.size (); ++j) {
				// 		elements [j]->recv_error ();
				// 	}
				// }
				for (int j = 0; j < (int) elements.size (); ++j) {
				// 	elements [j]->attempt_update ();
				// 	elements [j]->update ();
					elements [j]->update_timestep (t_timestep);
				}
			}
		}

	protected:
		int id; //!< The integer processor id
		int tsteps; //!< The integer total number of timesteps to take
		int count;
		
		std::vector <double> global_matrix;
		std::vector <double> global_rhs;
		std::vector <int> ipiv;
		std::vector <std::shared_ptr <element>> elements; //!< A vector containing shared pointers to the contained elements
		io::parameter_map inputParams; //!< The parameter map object containing the input parameters
		utils::messenger* messenger_ptr;
	};
} /* bases */
