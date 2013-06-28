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

namespace bases
{
	class element;
} /* bases */

namespace one_d
{
   /*!*******************************************************************
    * \brief An abstract plan class to handle 1D boundary conditions
    * 
    * This class implements boundary conditions in 1D with OpenMPI
    *********************************************************************/
	class mpi_boundary : public bases::boundary
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer width of the boundary
		 * \param i_ext_send The integer send tag
		 * \param i_ext_recv The integer recv tag
		 * \copydoc bases::boundary::boundary ()
		 *********************************************************************/
		mpi_boundary (bases::element* i_element_ptr, int i_edge, int i_ext_send, int i_ext_recv, int i_process) : bases::boundary (i_element_ptr, i_edge, i_process) {
			MTRACE ("Instantiating...");
			ext_send = i_ext_send;
			ext_recv = i_ext_recv;
			MTRACE ("Instantiated.");
		}
	
		virtual ~mpi_boundary () {}
		
		/*!*******************************************************************
		 * \copydoc bases::boundary::send ()
		 *********************************************************************/
		virtual void send (int name);
		
		/*!*******************************************************************
		 * \copydoc bases::boundary::recv ()
		 *********************************************************************/
		virtual void recv (int name);
				
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		virtual void execute () {
		}
		
	protected:
		int n; //!< The integer width of the boundary

		int ext_send; //!< The integer send tag
		int ext_recv; //!< The integer recv tag
		
		double buffer; //!< The buffer that contains the data
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
