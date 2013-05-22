/*!***********************************************************************
 * \file bases/boundary.hpp
 * Spectral Element
 * 
 * This file contains the base boundary class.
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

/*!*******************************************************************
 * \brief A set of flags to be used for boundaries
 * 
 * Currently, this only contains 1D boundaries, but more can easily be
 * added.
 *********************************************************************/
enum boundary_flags {
	linked_0 = 0x20,
	linked_n = 0x40
};

namespace bases
{
	class element;
	
	/*!*******************************************************************
	 * \brief A plan that implements a boundary condition
	 * 
	 * This plan adds the send and recv methods, which are used to
	 * communicate between boundaries with external memory. In general,
	 * this is done with a message passing interface, such as OpenMPI.
	 *********************************************************************/
	class boundary : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_edge An integer flag that specifies which boundary flag to use
		 * \param i_process The integer process id destination for send and recv
		 *********************************************************************/
		boundary (int i_edge, int i_process) : plan () {
			edge = i_edge;
			process = i_process;
		}
		
		virtual ~boundary () {}
		
		/*!*******************************************************************
		 * \copydoc plan::associate ()
		 *********************************************************************/
		virtual void associate (element* i_element_ptr);
		
		/*!*******************************************************************
		 * \brief Send the boundary information to the linked processor
		 *********************************************************************/
		virtual void send () {}
		
		/*!*******************************************************************
		 * \brief Receive the boundary information from the linked processor
		 *********************************************************************/
		virtual void recv () {}
	
	protected:
		int edge; //!< An integer flag that specifies which boundary to use
		int index; //!< The integer index to the boundary in the element
		int increment; //!< The integer increment to the next most inward index
		int process; //!< The integer process id destination for send and recv
	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_70XV2K58 */
