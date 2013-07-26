/*!***********************************************************************
 * \file bases/message.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef MESSENGER_HPP_JNOTV271
#define MESSENGER_HPP_JNOTV271

#include <vector>

namespace utils
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
		messenger (int* argc, char*** argv);
		
		virtual ~messenger ();
		
		virtual double& operator[] (int i);
		
		virtual void send (double* data, int process, int tag, double weight = 1.0, int size = 1, int inc = 1);
		
		virtual void recv (double* data, int process, int tag, double weight = 0.0, int size = 1, int inc = 1);
		
		virtual void min (double* data);
		
		virtual bool bool_and (bool boolean);
		
		int get_np () {
			return np;
		}
		
		int get_id () {
			return id;
		}
		
	private:
		int np;
		int id;
		
		std::vector <double> buffer;
		std::vector <int> int_buffer;
	};
} /* utils */

#endif /* end of include guard: MESSENGER_HPP_JNOTV271 */
