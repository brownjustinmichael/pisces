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
	 * \brief An implementation of the boundary class in 1D
	 * 
	 * This class handles the single point boundary by taking a weighted 
	 * average of the edges of two elements. It can also be used to zero 
	 * a boundary.
	 *********************************************************************/
	class boundary : public bases::boundary
	{
	public:
		/*!*******************************************************************
		 * \copydoc bases::boundary::boundary ()
		 *********************************************************************/
		boundary (int i_edge = fixed_0, bases::element* i_ext_element_ptr = NULL, int i_ext_edge = fixed_0) : bases::boundary (i_edge, i_ext_element_ptr, i_ext_edge) {}
	
		virtual ~boundary () {}
		
		virtual void fix_edge (int name) {
			fix_edge (name, (&((*element_ptr) [name])) [index]);
		}
		
		virtual void fix_edge (int name, double value) {
			fixed_points [name] = value;
		}
	
		/*!*******************************************************************
		 * \copydoc bases::boundary::execute ()
		 *********************************************************************/
		inline virtual void execute () {
			TRACE (logger, "Executing...");
			
			bases::boundary::execute ();
			
			for (std::map <int, double>::iterator iter = fixed_points.begin (); iter != fixed_points.end (); ++iter) {
				(&((*element_ptr) [iter->first])) [index] = iter->second;
				if (ext_element_ptr) {
					(&((*ext_element_ptr) [iter->first])) [ext_index] = iter->second;
				}
			}
			
			/*
				TODO Make implementation more general
			*/
		
			TRACE (logger, "executed.")
		}
	private:
		std::map <int, double> fixed_points;
	};
} /* one_d */


#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
