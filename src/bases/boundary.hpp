/*!***********************************************************************
 * \file bases/boundary.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_70XV2K58
#define BOUNDARY_HPP_70XV2K58

#include <vector>
#include <functional>
#include <map>
#include "../config.hpp"
#include "plan.hpp"

enum boundary_flags {
	fixed_0 = 0x20,
	fixed_n = 0x40
};

namespace bases
{
	class element;
	
   /*!*******************************************************************
    * \brief An abstract plan class to handle boundary conditions
    * 
    * This class contains the necessary coefficients and pointers to 
    * handle boundary conditions. The implementation will depend on the 
    * number of dimensions.
    *********************************************************************/
	class boundary : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_alpha_plus A double coefficient for the contribution from the positive boundary
		 * \param i_data_plus A pointer to the double first element of the positive boundary
		 * \param i_alpha_minus A double coefficient for the contribution from the negative boudary
		 * \param i_data_minus A pointer to the double first element of the negative boundary
		 *********************************************************************/
		boundary (int i_n, int i_to_next, int i_edge, bases::element* i_ext_element_ptr, int i_ext_edge) : plan () {
			MTRACE ("Instantiating...");
			n = 1;
			increase_n (i_n);
			to_next = i_to_next;
			edge = i_edge;
			ext_edge = i_ext_edge;
			ext_element_ptr = i_ext_element_ptr;
			MTRACE ("Instantiated.");
		}
	
		virtual ~boundary () {}
		
		virtual void associate (bases::element* i_element_ptr);
		
		virtual void increase_n (int i_n) {
			if (i_n > n) {
				n = i_n;
				send_buffer.resize (n);
				recv_buffer.resize (n);
			}
		}
		
		virtual void add_plan (plan* i_plan) {}
		
		virtual void set_ext_buffer (double* i_ext_buffer) {
			ext_buffer = i_ext_buffer;
		}
		
		virtual double* get_buffer () {
			return &send_buffer [0];
		}
		
		virtual void fix_edge (int name) = 0;
		
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () {}
		
		virtual void send ();
		
		virtual void recv ();
		
		friend class element;

	protected:
		int n;
		int to_next;
		int edge;
		int ext_edge;
		int index;
		int ext_index;
		double* ext_buffer;
		bases::element* ext_element_ptr;
		
		std::vector <double> send_buffer;
		std::vector <double> recv_buffer;
	};
	
	class active_boundary : public boundary
	{
	public:
		active_boundary (int i_n, int i_to_next, int i_edge, bases::element* i_ext_element_ptr, int i_ext_edge) : boundary (i_n, i_to_next, i_edge, i_ext_element_ptr, i_ext_edge) {
			n_plans = 0;
		}
		
		virtual ~active_boundary () {}
		
		virtual void add_plan (plan* i_plan) {
			if (i_plan) {
				plans.push_back (i_plan);
				++n_plans;
			}
		}
		
		virtual void execute ();
	
	protected:
		int n_plans;
		std::vector <plan*> plans;
	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_70XV2K58 */
