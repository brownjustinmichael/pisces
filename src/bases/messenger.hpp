/*!***********************************************************************
 * \file bases/messenger.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef MESSENGER_HPP_JNOTV271
#define MESSENGER_HPP_JNOTV271

#include <cassert>
#include <algorithm>
#include <vector>
#ifdef _MPI
#include "mpi.h"
#endif // _MPI
#include "../config.hpp"
#include "../utils/utils.hpp"

/*!**********************************************************************
 * An enum of various MPI actions that can be taken
 ************************************************************************/
enum modes {
	send_mode = 0,
	recv_mode = 1
};

namespace bases
{	
	/*!***********************************************************************
	 * \brief A class to manage boundary communication
	 * 
	 * This class should provide implementation to send boundary information
	 * to other elements. It may be dimension specific. There should only be 
	 * one per thread.
	 ************************************************************************/
	class messenger
	{
	public:
		/*!**********************************************************************
		 * \param argc An integer pointer to the number of command arguments
		 * \param argv A pointer to an array of character arrays of command arguments
		 * \param n_boundaries The integer number of boundaries per element (must be even)
		 ************************************************************************/
		messenger (int* argc, char*** argv, int n_boundaries);
		
		virtual ~messenger ();
		
		/*!**********************************************************************
		 * \brief Return the id of the current process
		 * 
		 * \return The integer id of the current process
		 ************************************************************************/
		int get_id () {
			return id;
		}
		
		/*!**********************************************************************
		 * \brief Return the total number of processes
		 * 
		 * \return The integer number of processes
		 ************************************************************************/
		int get_np () {
			return np;
		}
		
		/*!**********************************************************************
		 * \brief Return whether the specified edge is linked
		 * 
		 * \return A bool of if the edge is linked to another process
		 ************************************************************************/
		bool linked (int edge) {
			return process_queue [edge_to_index (send_mode, edge)] != -1;
		}
		
		/*!**********************************************************************
		 * \brief Given an edge and action mode, determine the correct index
		 * 
		 * Each action on each edge has a queue index, which this function reads
		 * 
		 * \param mode The integer representation of the action (0 for send, 1 for recv)
		 * \param edge The integer representation of the edge
		 ************************************************************************/
		virtual inline int edge_to_index (int mode, int edge) {
			if (id % 2 == 0) {
				if (edge % 2 == 0) {
					return 2 * edge + mode;
				} else {
					return 2 * edge + 1 - mode;
				}
			} else {
				if (edge % 2 == 0) {
					return 2 * (1 + edge) + mode;
				} else {
					return 2 * (edge - 1) + 1 - mode;
				}
			}
		}
	
		/*!**********************************************************************
		 * \brief Given an index, determine whether to send or recv
		 * 
		 * \param index The queue index of the action
		 ************************************************************************/
		virtual inline int index_to_mode (int index) {
			if (id % 2 == 0) {
				if (index % 4 == 0 || index % 4 == 3) {
					return send_mode;
				} else {
					return recv_mode;
				}
			} else {
				if (index % 4 == 1 || index % 4 == 2) {
					return send_mode;
				} else {
					return recv_mode;
				}
			}
		}

		/*!**********************************************************************
		 * \brief Link this element with another
		 * 
		 * \param edge The integer representation of the edge to link
		 * \param process The integer id of the process to link to
		 ************************************************************************/
		virtual void add_boundary (int edge, int process);
				
		/*!**********************************************************************
		 * \brief Add send action in the datatype queue
		 * 
		 * Warning. This does not send the data immediately. This could be 
		 * rewritten to save the information in a buffer, but it does not at 
		 * present. It sends if it can, but the information is not guaranteed to
		 * be sent until an entire set of send/recv commands are queued. 
		 * 
		 * \param n The integer number of elements to send
		 * \param data The datatype pointer to the data to send
		 * \param edge The integer representation of the corresponding edge
		 ************************************************************************/
		template <class datatype>
		void send (int n, datatype* data, int edge);
		
		template <class datatype>
		void send (int n, datatype* data, int process, int tag);
		
		/*!**********************************************************************
		 * \brief Add recv action in the datatype queue
		 * 
		 * Warning. This does not recv the data immediately. This could be 
		 * rewritten to wait for the information, but it does not at 
		 * present. It recvs if it can, but the information is not guaranteed to
		 * arrive until an entire set of send/recv commands are queued. 
		 * 
		 * \param n The integer number of elements to recv
		 * \param data The datatype pointer to the data to recv
		 * \param edge The integer representation of the corresponding edge
		 ************************************************************************/
		template <class datatype>
		void recv (int n, datatype* data, int edge);
		
		template <class datatype>
		void recv (int n, datatype* data, int process, int tag);
		
		template <class datatype>
		void gather (int n, datatype* data_in, datatype* data_out = NULL);
		
		template <class datatype>
		void gatherv (int n, datatype* data_in, int* ns, datatype* data_out = NULL);
		
		template <class datatype>
		void scatter (int n, datatype* data_in, datatype* data_out = NULL);
		
		template <class datatype>
		void scatterv (int n, datatype* data_in, int* ns, datatype* data_out = NULL);
		
		/*!**********************************************************************
		 * \brief Calculate a minimum across elements
		 * 
		 * This uses the messenger buffer.
		 * 
		 * \param data The datatype pointer to the datatype to minimize
		 ************************************************************************/	
		template <class datatype>	
		void min (datatype* data);
		
		/*!**********************************************************************
		 * \brief Determine if all processes meet a condition
		 * 
		 * This uses the messenger integer buffer.
		 * 
		 * \param boolean The bool condition to check
		 ************************************************************************/
		virtual bool bool_and (bool boolean);
		
	private:
		/*!**********************************************************************
		 * \brief Check to see if an action in the datatype queue is ready to be sent
		 ************************************************************************/
		template <class datatype>
		void data_check ();
		
		int id; //!< The integer id of the current process
		int np; //!< The integer number of total processes
		int data_iter; //!< The integer current location in the datatype queue
		int int_iter; //!< The integer current location in the integer queue
		
		std::vector <void*> data_queue; //!< A datatype pointer vector containing the elements to send/recv
		std::vector <int*> int_data_queue; //!< An integer pointer vector containing the elements to send/recv
		
		std::vector <int> n_queue; //!< An integer vector containing the number of elements to send/recv
		std::vector <int> process_queue; //!< An integer vector containing the processes for send/recv
	};
} /* bases */

#endif /* end of include guard: MESSENGER_HPP_JNOTV271 */
