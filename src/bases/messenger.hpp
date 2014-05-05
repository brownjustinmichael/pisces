/*!***********************************************************************
 * \file bases/messenger.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef MESSENGER_HPP_JNOTV271
#define MESSENGER_HPP_JNOTV271

#ifdef _MPI
#include "mpi.h"
#endif // _MPI
#include <cassert>
#include <algorithm>
#include <vector>
#include "../config.hpp"
#include "../utils/utils.hpp"

/*!**********************************************************************
 * An enum of various MPI actions that can be taken
 ************************************************************************/
enum modes {
	send_mode = 0,
	recv_mode = 1
};

enum mpi_flags {
	mpi_all_clear = 0x00,
	mpi_fatal = 0x01,
	mpi_skip = 0x02
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
		messenger (int* argc, char*** argv);
		
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
		
		void check_all (int *flags);

		void kill_all ();
		
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
		void recv (int n, datatype* data, int process, int tag);
		
		template <class datatype>
		void gather (int n, datatype* data_in, datatype* data_out = NULL);
		
		template <class datatype>
		void bcast (int n, datatype* data_in);
		
		template <class datatype>
		void gatherv (int n, datatype* data_in, int* ns, datatype* data_out = NULL);
		
		template <class datatype>
		void allgather (int n, datatype *data_in, datatype *data_out = NULL);
		
		template <class datatype>
		void scatter (int n, datatype* data_in, datatype* data_out = NULL);
		
		template <class datatype>
		void scatterv (int n, datatype* data_in, int* ns, datatype* data_out = NULL);
		
		template <class datatype>
		void allgatherv (int n, datatype* data_in, int* ns, datatype* data_out = NULL);
		
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
		
		std::vector <int> stati;
		
		int id; //!< The integer id of the current process
		int np; //!< The integer number of total processes
	};
} /* bases */

#endif /* end of include guard: MESSENGER_HPP_JNOTV271 */
