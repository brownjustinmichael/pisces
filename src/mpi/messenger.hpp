/*!***********************************************************************
 * \file utils/messenger.hpp
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
#include <typeinfo>

#include "versions/version.hpp"
#include "logger/logger.hpp"
#include "linalg/utils.hpp"

/*!**********************************************************************
 * An enum of MPI flags for checking the status of each process
 ************************************************************************/
enum mpi_flags {
	mpi_all_clear = 0x00,
	mpi_fatal = 0x01,
	mpi_skip = 0x02
};

/*!**********************************************************************
 * \namespace mpi
 * 
 * This namespace contains classes and functions relevant to MPI communication, such as the messenger class.
 ************************************************************************/
namespace mpi
{	
#ifdef _MPI
	/*!**********************************************************************
	 * \brief A helper function for MPI typing, for convenience
	 * 
	 * \param type The type_info object describing the desired MPI type
	 * 
	 * \return The MPI type associated with the argument
	 ************************************************************************/
	inline MPI_Datatype mpi_type (const std::type_info* type) {
		if (type == &typeid (double)) {
			return MPI::DOUBLE;
		} else if (type == &typeid (int)) {
			return MPI::INT;
		} else if (type == &typeid (float)) {
			return MPI::REAL;
		} else {
			FATAL ("Unrecognized MPI type");
			throw 0;
		}
	}
	
#endif // _MPI
	
	/*!***********************************************************************
	 * \brief A class to manage boundary communication
	 * 
	 * This class should provide implementation to send boundary information to other elements. It may be dimension specific. There should only be one per thread.
	 ************************************************************************/
	class messenger
	{
	private:
		int id; //!< The integer id of the current process
		int np; //!< The integer number of total processes
		
		std::vector <int> stati; //!< A vector of the status of each processor
		
		class config
		{
		public:
			config (int *argc = NULL, char *** argv = NULL) {
				// Initialize MPI
				MPI::Init (*argc, *argv);
			}
			
			virtual ~config () {
				// Finalize MPI
				MPI::Finalize ();
			}
		};
		
	public:
		/*!**********************************************************************
		 * \param argc An integer pointer to the number of command arguments
		 * \param argv A pointer to an array of character arrays of command arguments
		 ************************************************************************/
		messenger (int* argc = NULL, char*** argv = NULL) {
#ifdef _MPI
			static config config_instance (argc, argv);
			np = MPI::COMM_WORLD.Get_size ();
			id = MPI::COMM_WORLD.Get_rank ();
			
			// Enable the MPI error handler
			MPI::COMM_WORLD.Set_errhandler(MPI::ERRORS_THROW_EXCEPTIONS);
#else
			np = 1;
			id = 0;
#endif // _MPI
			stati.resize (np);
		}
		
		virtual ~messenger () {
			kill_all ();
		}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.0.2.0");
			return version;
		}
		
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
		 * \brief Check that all processors are working
		 * 
		 * This method checks that all MPI processes are behaving normally. If there's a fatal issue with MPI, this method throws an error.
		 * 
		 * \return True if all clear, false if the next process should be skipped
		 ************************************************************************/
		bool check_all () {
			TRACE ("Checking...");
			int flags = mpi_all_clear;
#ifdef _MPI
			// Gather the status of each messenger across the processes
			MPI::COMM_WORLD.Gather (&flags, 1, mpi_type (&typeid (int)), &stati [0], 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
			// Combine all the statuses into one set of binary flags
			/*
				TODO Because this is small, it's probably cheaper to combine the flags on every processor
			*/
			if (id == 0) {
				for (int i = 0; i < np; ++i) {
					flags |= stati [i];
				}
			}
#ifdef _MPI
			// Broadcast the combined flags
			MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
			// If there's a fatal flag, kill all processes
			if (flags & mpi_fatal) {
				throw 0;
			}
			TRACE ("Check complete.");
			// Return whether every processor is ready for instructions
			if (flags != mpi_all_clear) {
				return false;
			} else {
				return true;
			}
		}
		
		/*!**********************************************************************
		 * \brief Force the computer to skip any active MPI commands
		 ************************************************************************/
		void skip_all () {
			int flags = mpi_skip;
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
		}
		
		/*!**********************************************************************
		 * \brief Force the computer to kill all MPI processes
		 ************************************************************************/
		void kill_all () {
			int flags = mpi_fatal;
#ifdef _MPI
			MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
		}
		
		/*!**********************************************************************
		 * \brief Send the data to a given process with a specific tag
		 * 
		 * \param n The integer number of elements to send
		 * \param data The datatype pointer to the data to send
		 * \param process The integer rank of the receiving processor
		 * \param tag The identifying tag on the message
		 ************************************************************************/
		template <class datatype>
		void send (int n, const datatype* data, int process, int tag) {
#ifdef _MPI
			MPI::COMM_WORLD.Send (data, n, mpi_type (&typeid (datatype)), process, tag);
#endif // _MPI
		}
		
		/*!**********************************************************************
		 * \brief Receive the data from a given process with a specific tag
		 * 
		 * \param n The integer number of elements to recieve
		 * \param data The datatype pointer to the data to receive
		 * \param process The integer rank of the sending processor
		 * \param tag The identifying tag on the message
		 ************************************************************************/
		template <class datatype>
		void recv (int n, datatype* data, int process, int tag) {
		/*
			TODO Make check_two
		*/
#ifdef _MPI
			MPI::COMM_WORLD.Recv (data, n, mpi_type (&typeid (datatype)), process, tag);
#endif // _MPI
		}
		
		/*!**********************************************************************
		 * \brief Gather arrays from all processes onto process 0
		 * 
		 * \param n The number of elements to gather per process
		 * \param data_in The array of data to be gathered
		 * \param data_out The output array of all the data on process 0 (can be the same as data_in)
		 ************************************************************************/
		template <class datatype>
		void gather (int n, const datatype* data_in, datatype* data_out) {
			if (check_all ()) {
#ifdef _MPI
				if (id == 0 && data_out == data_in) {
					// If the gathering is in place, you need a spectial flag
					MPI::COMM_WORLD.Gather (MPI_IN_PLACE, n, mpi_type (&typeid (datatype)), data_out, n, mpi_type (&typeid (datatype)), 0);
				} else {
					MPI::COMM_WORLD.Gather (data_in, n, mpi_type (&typeid (datatype)), data_out, n, mpi_type (&typeid (datatype)), 0);
				}
#endif // _MPI
			}
		}
		
		/*!**********************************************************************
		 * \brief Gather arrays of varying length from all processes onto process 0
		 * 
		 * \param n The number of elements to send from an individual process
		 * \param data_in The input array of data to be gathered
		 * \param ns The array of the lengths of the arrays across all processes
		 * \param data_out The output array of all the data on process 0 (can be the same as data_in)
		 ************************************************************************/
		template <class datatype>
		void gatherv (int n, const datatype* data_in, int* ns, datatype* data_out) {
			if (check_all ()) {
#ifdef _MPI
				// Create a displacement vector that tracks the starting position of each gathered array
				std::vector <int> displs;
				if (id == 0) {
					displs.resize (np);
					for (int i = 1; i < np; ++i) {
						displs [i] = displs [i - 1] + ns [i - 1];
					}
				}
				if (id == 0 && data_out == data_in) {
					// If the gathering is in place, you need a spectial flag
					MPI::COMM_WORLD.Gatherv (MPI_IN_PLACE, n, mpi_type (&typeid (datatype)), data_in, ns, &displs [0], mpi_type (&typeid (datatype)), 0);
				} else {
					MPI::COMM_WORLD.Gatherv (data_in, n, mpi_type (&typeid (datatype)), data_out, ns, &displs [0], mpi_type (&typeid (datatype)), 0);
				}
#endif // _MPI
			}
		}
		
		/*!**********************************************************************
		 * \brief Gather arrays from all processes onto all processes
		 * 
		 * \param n The number of elements to gather per process
		 * \param data_in The array of data to be gathered
		 * \param data_out The output array of all the data (can be the same as data_in)
		 ************************************************************************/
		template <class datatype>
		void allgather (int n, const datatype *data_in, datatype *data_out) {
			if (check_all ()) {
#ifdef _MPI
				if (data_out == data_in) {
					// If the gathering is in place, you need a spectial flag
					MPI::COMM_WORLD.Allgather (MPI_IN_PLACE, n, mpi_type (&typeid (datatype)), data_out, n, mpi_type (&typeid (datatype)));
				} else {
					MPI::COMM_WORLD.Allgather (data_in, n, mpi_type (&typeid (datatype)), data_out, n, mpi_type (&typeid (datatype)));
				}
#endif // _MPI
			}
		}
		
		/*!**********************************************************************
		 * \brief Gather arrays of variable length from all processes onto all processes
		 * 
		 * \param n The number of elements to gather for each process
		 * \param data_in The array of data to be gathered
		 * \param ns The integer array of the lengths of the data on each processor
		 * \param data_out The output array of all the data (can be the same as data_in)
		 ************************************************************************/
		template <class datatype>
		void allgatherv (int n, const datatype* data_in, int* ns, datatype* data_out) {
			if (check_all ()) {
#ifdef _MPI
				// Create a displacement vector that tracks the starting position of each gathered array
				std::vector <int> displs;
				displs.resize (np);
				for (int i = 1; i < np; ++i) {
					displs [i] = displs [i - 1] + ns [i - 1];
				}
				if (data_out == data_in) {
					// If the gathering is in place, you need a spectial flag
					MPI::COMM_WORLD.Allgatherv (MPI_IN_PLACE, n, mpi_type (&typeid (datatype)), data_out, ns, &displs [0], mpi_type (&typeid (datatype)));
				} else {
					MPI::COMM_WORLD.Allgatherv ((const void *) data_in, n, mpi_type (&typeid (datatype)), (void *) data_out, ns, &displs [0], mpi_type (&typeid (datatype)));
				}
#endif // _MPI
			}
		}
		
		/*!**********************************************************************
		 * \brief Broadcast array from process 0 to all processes
		 * 
		 * \param n The number of elements to broadcast
		 * \param data The array of data to be broadcast on input and the final output array
		 ************************************************************************/
		template <class datatype>
		void bcast (int n, datatype* data) {
			if (check_all ()) {
#ifdef _MPI
				MPI::COMM_WORLD.Bcast (data, n, mpi_type (&typeid (datatype)), 0);
#endif // _MPI
			}
		}
		
		/*!**********************************************************************
		 * \brief Calculate a minimum across elements
		 * 
		 * This uses a messenger buffer.
		 * 
		 * \param data The datatype pointer to the datatype to minimize
		 ************************************************************************/	
		template <class datatype>	
		void min (datatype* data) {
			if (check_all ()) {
				if (np != 1) {
					std::vector <datatype> buffer (np);
#ifdef _MPI
					MPI::COMM_WORLD.Gather (data, 1, mpi_type (&typeid (datatype)), &buffer [0], 1, mpi_type (&typeid (datatype)), 0);
#endif // _MPI
					if (id == 0) {
						*data = *std::min_element (buffer.begin (), buffer.end ());
					}
#ifdef _MPI
					MPI::COMM_WORLD.Bcast (data, 1, mpi_type (&typeid (datatype)), 0);
#endif // _MPI
				}
			}
		}
	};
} /* mpi */

#endif /* end of include guard: MESSENGER_HPP_JNOTV271 */
