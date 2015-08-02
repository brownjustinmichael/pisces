/*!***********************************************************************
 * \file messenger.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/
 
 #include "messenger.hpp"
 
 namespace mpi
 {	
 	class config
 	{
	protected:
		bool finalize;

 	public:
 		config (int *argc = NULL, char *** argv = NULL, bool i_finalize = true) : finalize (i_finalize) {
 			// Initialize MPI
 			#ifdef _MPI
 			MPI::Init (*argc, *argv);
 			#endif
 		}
 		
 		virtual ~config () {
 			// Finalize MPI
 			#ifdef _MPI
 			if (finalize) MPI::Finalize ();
 			#endif
 		}
 	};

	messenger::messenger (int* argc, char*** argv, bool finalize) {
#ifdef _MPI
		static config config_instance (argc, argv, finalize);
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
	
	bool messenger::check_all () {
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
	
	void messenger::skip_all () {
#ifdef _MPI
			int flags = mpi_skip;
			MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
	}
	
	void messenger::kill_all () {
#ifdef _MPI
			int flags = mpi_fatal;
			MPI::COMM_WORLD.Bcast (&flags, 1, mpi_type (&typeid (int)), 0);
#endif // _MPI
	}
} /* mpi */
