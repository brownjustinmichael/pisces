/*!***********************************************************************
 * \file one_d/element.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_3SURDTOH
#define ELEMENT_HPP_3SURDTOH

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "../config.hpp"
#include "../bases/element.hpp"
#include "../bases/plan.hpp"
	
namespace one_d
{
	/*!*******************************************************************
	 * A 1D implementation of the element base class. This provides the
	 * storage and indexing facilities. The plans should be added in a 
	 * further subclass.
	 * 
	 * \brief \copybrief bases::element
	 *********************************************************************/
	class element : public bases::element
	{
	public:	
		/*!*******************************************************************
		 * \param i_n The number of data elements in each scalar
		 * \copydoc bases::element::element ()
		 *********************************************************************/
		element (int i_n, int i_flags) : bases::element (i_flags) {
			n = i_n;
			
			cell.resize (i_n);
			for (int i = 0; i < i_n; ++i) {
				cell [i] = i;
			}
		}
		
		virtual ~element () {}
		
		/*!*******************************************************************
		 * \copydoc bases::element::add_scalar ()
		 *********************************************************************/
		inline void add_scalar (int name) {
			scalars [name].resize (n + 1);
		}
		
		inline void execute_boundaries () {
			// for (int i = 0; i < n; ++i) {
			// 	DEBUG (logger, "rhs [" << i << "] = " << (*this) (rhs, i));
			// }

			bases::element::execute_boundaries ();
			
			// for (int i = 0; i < n; ++i) {
			// 	DEBUG (logger, "rhs [" << i << "] = " << (*this) (rhs, i));
			// }
		}
		
		
		/*!*******************************************************************
		 * \copydoc bases::element::operator[] ()
		 *********************************************************************/
		inline double &operator[] (int name) {
			return scalars [name] [0];
		}
		
	protected:
		int n; //!< The number of elements in each 1D array
	   	std::vector<int> cell; //!< An integer array for tracking each cell number for output

		std::map<int, std::vector<double>> scalars; //!< A vector of scalar vectors
	};

	/*!*******************************************************************
	 * \brief A simple implementation of the element class with diffusion
	 * 
	 * This class contains a full element's capacity to run a single 
	 * element diffusion in 1D. It cannot yet tie to other elements or 
	 * calculate its own timestep. 
	 *********************************************************************/
	class advection_diffusion_element : public element
	{
	public:
		/*!*******************************************************************
		 * \param i_n The number of elements in each 1D data array
		 * \param i_flags Flags for the boundary conditions and evaluation
		 *********************************************************************/
		advection_diffusion_element (std::string name, double i_alpha_plus, double i_alpha_minus, int i_n, double *initial_position, double *intial_velocity, int i_flags);
		virtual ~advection_diffusion_element () {}
		
	private:
		double alpha_plus;
		double alpha_minus;
		std::vector<double> matrix; //!< A vector containing the double matrix used in the implicit solver
	};
} /* one_d */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
