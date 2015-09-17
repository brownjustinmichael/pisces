/*!**********************************************************************
 * \file rezone.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-02-25.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REZONE_CPP_UEW78T7Q
#define REZONE_CPP_UEW78T7Q

#include "io/output.hpp"
#include "io/formats/virtual.hpp"
#include "linalg/interpolate.hpp"
#include "plans/grids/grid.hpp"
#include "linalg/utils.hpp"

namespace pisces
{
	template <class datatype>
	class element;
	
	/**
	 * @brief A collection of data and methods useful in rezoning data
	 * @details The rezone data class is designed to minimize the necessary interface with the gsl library and with the virtual file objects needed to rezone the data.
	 */
	template <class datatype>
	struct rezone_data {
		element <datatype> *element_ptr; //!< A pointer to the element to rezone
		int id; //!< The id of the process
		int np; //!< The number of processes currently running
		/*
			TODO It might be possible to make this dynamic
		*/
		datatype positions [64]; //!< The array of positions used in the iterations
		datatype min_size; //!< The minimmum extent of any element
		datatype max_size; //!< The maximum extent of any element
		double (*merit_func) (element <datatype> *element_ptr, formats::virtual_file *); //!< A merit function to be maximized
		
		/**
		 * @brief Given the merit function associated with the rezone data, this function is in the correct form for gsl to parse
		 * 
		 * @param i_rezone_data A pointer to the rezone data
		 * @return The merit value that the merit function has returned, summed over all elements
		 */
		static double func (void *i_rezone_data) {
			element <datatype> *element_ptr = ((rezone_data <datatype> *) i_rezone_data)->element_ptr;
			mpi::messenger *messenger_ptr = element_ptr->messenger_ptr;
			double value = ((rezone_data <datatype> *)i_rezone_data)->merit_func (element_ptr, element_ptr->make_rezoned_virtual_file (((rezone_data <datatype> *)i_rezone_data)->positions, &*(element_ptr->rezone_virtual_file), profile_only));
			messenger_ptr->sum (&value);
			return value;
		}
		
		/*!**********************************************************************
		 * \brief Given two sets of zoning data, calculate the distance in parameter space
		 * 
		 * \param i_new_rezone_data A pointer to an array of new rezone_union objects
		 * \param i_old_rezone_data A pointer to an array of old rezone_union objects
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static double rezone_step_size (void *i_new_rezone_data, void *i_old_rezone_data) {
			datatype total = 0;
			datatype *new_positions = ((rezone_data <datatype> *) i_new_rezone_data)->positions;
			datatype *old_positions = ((rezone_data <datatype> *) i_old_rezone_data)->positions;
			for (int i = 1; i < ((rezone_data <datatype> *) i_new_rezone_data)->np; ++i) {
				total += (new_positions [i] - old_positions [i]) * (new_positions [i] - old_positions [i]);
			}
			return sqrt (total);
		}

		/*!**********************************************************************
		 * \brief Generate a new set of parameters from the old ones with a random number generator
		 * 
		 * \param r A pointer to a gnu scientific library random number generator
		 * \param i_rezone_data A pointer to an array of rezone_union objects
		 * \param step_size The datatype step_size in parameter space
		 * 
		 * In order, the rezone_union objects must contain 1) a pointer to the element, 2) the total number of elements in the system, 3) min_size argument from rezone_minimize_ts, 4) max_size argument from rezone_minimize_ts, 5-n) positions of the zonal boundaries. This method is to be called by the GSL simulated annealing routine.
		 ************************************************************************/
		static void rezone_generate_step (const gsl_rng *r, void *i_rezone_data, datatype step_size) {
			element <datatype> *element_ptr = ((rezone_data <datatype> *) i_rezone_data)->element_ptr;
			datatype *old_positions = ((rezone_data <datatype> *) i_rezone_data)->positions;
			datatype min_size = ((rezone_data <datatype> *) i_rezone_data)->min_size;
			datatype max_size = ((rezone_data <datatype> *) i_rezone_data)->max_size;
			datatype positions [64];
			
			int id = ((rezone_data <datatype> *) i_rezone_data)->id;
			int np = ((rezone_data <datatype> *) i_rezone_data)->np;
			mpi::messenger *messenger_ptr = element_ptr->messenger_ptr;
			if (id == 0) {
				// Generate a random radius that is less than or equal to the step size
				datatype radius = gsl_rng_uniform (r) * step_size;
				if (radius == 0.0) {
					messenger_ptr->skip_all ();
					return;
				}
				// Generate a random step for every zonal position and sum the total step size
				datatype xs [np + 1];
				datatype total = 0.0;
				for (int i = 1; i < np; ++i) {
					xs [i] = gsl_rng_uniform (r) * 2.0 - 1.0;
					total += xs [i] * xs [i];
				}
				// If possible, rescale the total step to be equal to the radius; otherwise, change nothing
				if (total == 0.0) {
					messenger_ptr->skip_all ();
					return;
				}
				total = sqrt (total);
				total /= radius;
				positions [0] = old_positions [0];
				for (int i = 1; i < np; ++i) {
					positions [i] = xs [i] / total + old_positions [i];
				}
				positions [np] = old_positions [np];
				// Broadcast the new positions
				messenger_ptr->template bcast <datatype> (np + 1, positions);
			} else {
				// Receive the new positions
				messenger_ptr->template bcast <datatype> (np + 1, positions);
			}
			// Iterate forward through the data, checking that the mininum and maximum sizes are obeyed
			old_positions [0] = positions [0];
			for (int i = 1; i < np; ++i) {
				// Check minimum size
				if (positions [i] < old_positions [i - 1] + min_size) {
					positions [i] = old_positions [i - 1] + min_size;
				}
				// Check maximum size
				if (positions [i] > old_positions [i - 1] + max_size) {
					positions [i] = old_positions [i - 1] + max_size;
				}
				old_positions [i] = positions [i];
			}
			// Iterate backward through the positions, checking that the minimum and maximum sizes are obeyed
			for (int i = np - 1; i >= 1; --i) {
				// Check minimum size
				if (old_positions [i] > old_positions [i + 1] - min_size) {
					old_positions [i] = old_positions [i + 1] - min_size;
				}
				// Check maximum size
				if (old_positions [i] < old_positions [i + 1] - max_size) {
					old_positions [i] = old_positions [i + 1] - max_size;
				}
				DEBUG ("POS: " << old_positions [i]);
			}
			// Make sure we didn't accidentally run off the end of the position array
			if (old_positions [0] < old_positions [1] - max_size) {
				for (int i = 0; i < np; ++i) {
					FATAL (old_positions [i] << " " << max_size);
				}
				FATAL ("Unrealistic size constraints in rezone: max_size too small");
				throw 0;
			}
			if (old_positions [0] > old_positions [1] - min_size) {
				for (int i = 0; i < np; ++i) {
					FATAL (old_positions [i] << " " << min_size);
				}
				FATAL ("Unrealistic size constraints in rezone: min_size too large");
				throw 0;
			}
		}
	};
	
	/*!**********************************************************************
	 * \brief Rezone the elements according to the output grid
	 * 
	 * \param inter_messenger The mpi messenger that communicates between the elements
	 * \param input_grid A pointer to the grid describing the extent of the input data
	 * \param output_grid A pointer to the grid describing the desired extent of the output data
	 * \param input_virtual_file A pointer to the virtual file containing the input data
	 * \param output_virtual_file A pointer to the virtual file constructed by rezone. If NULL, use the input_virtual_file
	 * @param value_buffer A pointer to a buffer to hold the entire data of one variable
	 * @param inter_buffer A pointer to a buffer to hold the positional data of the vertical dimension
	 * 
	 * This method takes an input_virtual_file and input_grid and rezones them according to the extent of output_grid. It will do so by communicating with the other elements to collect the data necessary for the rezone.
	 ************************************************************************/
	template <class datatype>
	void rezone (mpi::messenger *inter_messenger, grids::grid <datatype> *input_grid, grids::grid <datatype> *output_grid, formats::virtual_file *input_virtual_file, formats::virtual_file *output_virtual_file, datatype *value_buffer, datatype *inter_buffer) {
		if (output_virtual_file != input_virtual_file) {
			*output_virtual_file = *input_virtual_file;
		}
		
		int nn = input_grid->get_n ();
		std::vector <int> ns (inter_messenger->get_np (), nn);
		std::vector <int> nns (inter_messenger->get_np (), nn);
		inter_messenger->allgather <int> (1, &ns [0], &ns [0]);
		int nsum = 0;
		int nhere = 0;
		for (int i = 0; i < inter_messenger->get_np (); ++i) {
			nns [i] = ns [i];
			nsum += ns [i];
		}
		std::vector <datatype> position_buffer (nsum, 0.0);
		
		// Gather the global positions
		/*
			TODO Move this outside of the loop
		*/
		inter_messenger->allgatherv <datatype> (nn, &((*input_grid) [0]), &ns [0], &position_buffer [0]);
		
		// Iterate through the data
		for (typename std::map <std::string, void *>::iterator iter = input_virtual_file->begin (); iter != input_virtual_file->end (); iter++) {
			if (input_virtual_file->dims [iter->first] [1] != 1 && input_virtual_file->check_type <datatype> (iter->first)) {
				TRACE ("Rezoning " << iter->first << "...");
				
				nhere = 0;
				for (int i = 0; i < inter_messenger->get_np (); ++i) {
					nns [i] = ns [i] * input_virtual_file->dims [iter->first] [0];
					if (i < inter_messenger->get_id ()) {
						nhere += nns [i];
					}
				}
				
				// Gather the entirety of the values for the data in every element; need to switch to row-major order first
				linalg::matrix_switch (nn, input_virtual_file->dims [iter->first] [0], &(input_virtual_file->index <datatype> (iter->first)), &inter_buffer [nhere]);
				inter_messenger->allgatherv <datatype> (nn * input_virtual_file->dims [iter->first] [0], inter_buffer, &nns [0], inter_buffer);
				linalg::matrix_switch (input_virtual_file->dims [iter->first] [0], nsum, inter_buffer, value_buffer);
				output_virtual_file->add_var <datatype> (iter->first, input_virtual_file->dims [iter->first] [0], input_grid->get_n ());
				
				// Interpolate the new values from the global positions and values
				linalg::interpolate <datatype> (output_grid->get_n (), output_virtual_file->dims [iter->first] [0], nsum, 1.0, 0.0, &position_buffer [0], value_buffer, &((*output_grid) [0]), &(output_virtual_file->index <datatype> (iter->first)));
			}
		}
	}
} /* pisces */


#endif /* end of include guard: REZONE_CPP_UEW78T7Q */
