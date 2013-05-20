/*!***********************************************************************
 * \file one_d/boundary.hpp
 * Spectral Element
 * 
 * This file provides the implementation for the boundary base class. 
 * New realizations of boundary conditions should include this file to 
 * inherit from the boundary class.
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_CP5DR4CP
#define BOUNDARY_HPP_CP5DR4CP

#include <memory>
#include <map>
#include "../bases/boundary.hpp"
#include "../bases/element.hpp"
#include "../config.hpp"

namespace one_d
{
   /*!*******************************************************************
    * \brief An abstract plan class to handle boundary conditions
    * 
    * This class contains the necessary coefficients and pointers to 
    * handle boundary conditions. The implementation will depend on the 
    * number of dimensions.
    *********************************************************************/
	class link_boundary : public bases::boundary
	{
	public:
		typedef std::map <int, std::vector <double>> buffer;
		
		/*!*******************************************************************
		 * \param i_alpha_plus A double coefficient for the contribution from the positive boundary
		 * \param i_data_plus A pointer to the double first element of the positive boundary
		 * \param i_alpha_minus A double coefficient for the contribution from the negative boudary
		 * \param i_data_minus A pointer to the double first element of the negative boundary
		 *********************************************************************/
		link_boundary (int i_n, int i_edge, double* i_ext_send, double* i_ext_recv) : bases::boundary (i_edge) {
			MTRACE ("Instantiating...");
			n = i_n;
			ext_send = i_ext_send;
			ext_recv = i_ext_recv;
			MTRACE ("Instantiated.");
		}
	
		virtual ~link_boundary () {}
		
		virtual void associate (bases::element* element_ptr) {
			bases::boundary::associate (element_ptr);
			for (bases::element::iterator iter = element_ptr->begin (); iter != element_ptr->end (); ++iter) {
				send_buffer [*iter].resize (n);
				recv_buffer [*iter].resize (n);
			}
		}
		
		virtual void send ();
		
		virtual void recv ();
				
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute () {
			bases::boundary::execute ();
			for (bases::element::iterator iter = element_ptr->begin (); iter != element_ptr->end (); ++iter) {
				MDEBUG ("" << *iter << " send: " << send_buffer [*iter] [0] << " recv: " << recv_buffer [*iter] [0]);
				(*element_ptr) (*iter, index) = 0.5 * send_buffer [*iter] [0] + 0.5 * recv_buffer [*iter] [0];
			}	
		}
		
	protected:
		int n;

		double* ext_send;
		double* ext_recv;
		
		buffer send_buffer;
		buffer recv_buffer;
	};
	
	/*!*******************************************************************
	 * \brief An implementation of the boundary class in 1D
	 * 
	 * This class handles the single point boundary by taking a weighted 
	 * average of the edges of two elements. It can also be used to zero 
	 * a boundary.
	 *********************************************************************/
	class diffusive_boundary : public link_boundary
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::link_boundary::link_boundary ()
		 *********************************************************************/
		diffusive_boundary (int i_edge, double i_coeff, int i_position, int i_name_in, int i_name_out, double* i_ext_send, double* i_ext_recv) : link_boundary (2, i_edge, i_ext_send, i_ext_recv) {
			coeff = i_coeff;
			position = i_position;
			name_in = i_name_in;
			name_out = i_name_out;
		}
	
		virtual ~diffusive_boundary () {}
		
		/*!*******************************************************************
		 * \copydoc bases::link_boundary::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE (logger, "Executing...");
			double x01 = (send_buffer [position] [1] - send_buffer [position] [0]);
			double x02 = (send_buffer [position] [1] - recv_buffer [position] [1]);
			double x12 = (send_buffer [position] [0] - recv_buffer [position] [1]);
			
			send_buffer [name_out] [0] += 2.0 * coeff * 2.0 * (send_buffer [name_in] [1] * x12 - send_buffer [name_in] [0] * x02 + recv_buffer [velocity] [1] * x01) / x01 / x02 / x12;
			recv_buffer [name_out] [0] += send_buffer [name_in] [0];
			
			link_boundary::execute ();
			TRACE (logger, "executed.")
		}
		
	protected:
		int position;
		int name_in;
		int name_out;
		double coeff;
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
