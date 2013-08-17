/*!***********************************************************************
 * \file bases/message.hpp
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
#include "mpi.h"
#include "../config.hpp"
#include "../utils/utils.hpp"

enum mode {
	send_mode = 0,
	recv_mode = 1
};

namespace bases
{	
	/*!***********************************************************************
	 * This class should provide implementation to send boundary information
	 * to other elements. It may be dimension specific. As there is only need
	 * for one of these per element, it may be natural to incorporate it, but
	 * that is not certain.
	 ************************************************************************/
	class messenger
	{
	public:
		messenger (int* argc, char*** argv, int n_boundaries);
		
		virtual ~messenger ();
		
		virtual void add_boundary (int edge, int process);
		
		virtual int edge_to_index (int mode, int edge);
		
		virtual int index_to_mode (int index);
		
		virtual void send (int n, double* data, int edge, int inc = 1);
		
		virtual void recv (int n, double* data, int edge, int inc = 1);

		virtual void double_check ();

		virtual void send (int n, int* data, int edge, int inc = 1);
		
		virtual void recv (int n, int* data, int edge, int inc = 1);
		
		virtual void int_check ();
		
		virtual void min (double* data);
		
		virtual bool bool_and (bool boolean);
		
		int get_np () {
			return np;
		}
		
		int get_id () {
			return id;
		}
		
		bool linked (int edge) {
			return process_queue [edge_to_index (send_mode, edge)] != -1;
		}
		
	private:
		int np;
		int id;
		int double_iter;
		int int_iter;
		
		std::vector <int> n_queue;
		std::vector <double*> double_data_queue;
		std::vector <int*> int_data_queue;
		std::vector <int> process_queue;
		
		std::vector <double> buffer;
		std::vector <int> int_buffer;
	};
} /* bases */

#endif /* end of include guard: MESSENGER_HPP_JNOTV271 */
