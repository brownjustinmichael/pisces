/*!***********************************************************************
 * \file boundary_one_d.hpp
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
    * \brief An abstract plan class to handle 1D boundary conditions
    * 
    * This class implements boundary conditions in 1D with OpenMPI
    *********************************************************************/
	class link_boundary : public bases::boundary
	{
	public:
		typedef std::map <int, std::vector <double>> buffer;
		
		/*!*******************************************************************
		 * \param i_n The integer width of the boundary
		 * \param i_ext_send The integer send tag
		 * \param i_ext_recv The integer recv tag
		 * \copydoc bases::boundary::boundary ()
		 *********************************************************************/
		link_boundary (int i_n, int i_edge, int i_ext_send, int i_ext_recv, int i_process) : bases::boundary (i_edge, i_process) {
			MTRACE ("Instantiating...");
			n = i_n;
			ext_send = i_ext_send;
			ext_recv = i_ext_recv;
			MTRACE ("Instantiated.");
		}
	
		virtual ~link_boundary () {}
		
		/*!*******************************************************************
		 * \copydoc bases::boundary::associate ()
		 *********************************************************************/
		virtual void associate (bases::element* element_ptr) {
			bases::boundary::associate (element_ptr);
			int j = 0;
			for (bases::element::iterator iter = element_ptr->begin (); iter != element_ptr->end (); ++iter) {
				send_buffer [*iter].resize (n);
				recv_buffer [*iter].resize (n);
				++j;
			}
			to_send.resize (n * j);
			to_recv.resize (n * j);
		}
		
		/*!*******************************************************************
		 * \copydoc bases::boundary::send ()
		 *********************************************************************/
		virtual void send ();
		
		/*!*******************************************************************
		 * \copydoc bases::boundary::recv ()
		 *********************************************************************/
		virtual void recv ();
				
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		virtual void execute () {
			bases::boundary::execute ();
			for (bases::element::iterator iter = element_ptr->begin (); iter != element_ptr->end (); ++iter) {
				(*element_ptr) (*iter, index) = 0.5 * send_buffer [*iter] [0] + 0.5 * recv_buffer [*iter] [0];
			}	
		}
		
	protected:
		int n; //!< The integer width of the boundary

		int ext_send; //!< The integer send tag
		int ext_recv; //!< The integer recv tag
		
		std::vector <double> to_send; //!< The serialized vector containing the data to send
		std::vector <double> to_recv; //!< The serialized vector to contain the received data
		
		buffer send_buffer; //!< The buffer that contains the data to send
		buffer recv_buffer; //!< The buffer that contains the received data
		
		/*
			TODO We have a redundancy here that should be fixed.
		*/
	};
	
	/*!*******************************************************************
	 * \brief An implementation of the boundary class in 1D for boundaries with viscosity
	 * 
	 * This class handles the single point boundary by using Lagrange 
	 * interpolation using two points from each element.
	 *********************************************************************/
	class diffusive_boundary : public link_boundary
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::link_boundary::link_boundary ()
		 *********************************************************************/
		diffusive_boundary (int i_edge, int i_ext_send, int i_ext_recv, int i_process) : link_boundary (2, i_edge, i_ext_send, i_ext_recv, i_process) {}
	
		/*!*******************************************************************
		 * \copydoc bases::link_boundary::associate ()
		 *********************************************************************/
		void associate (bases::element* i_element_ptr) {
			link_boundary::associate (i_element_ptr);
			coeff = element_ptr->get_dparam ("diffusion_coeff");
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
			
			send_buffer [rhs] [0] += 2.0 * coeff * 2.0 * (send_buffer [velocity] [1] * x12 - send_buffer [velocity] [0] * x02 + recv_buffer [velocity] [1] * x01) / x01 / x02 / x12;
			recv_buffer [rhs] [0] += send_buffer [rhs] [0];
			
			link_boundary::execute ();
			TRACE (logger, "executed.")
		}
		
	protected:
		double coeff; //!< The double diffusion coefficient
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
