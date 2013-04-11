/*!***********************************************************************
 * \file element.hpp
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
#include <fftw3.h>
#include "../diffusion/diffusion.hpp"
#include "../io/io.hpp"

namespace element
{
	enum index {
		position = 0,
		velocity = 1,
		rhs = 2
	};
	
	/*!*******************************************************************
	 * \brief This is the basic class of the code
	 * 
	 * A true run will contain multiple elements tied together at the 
	 * boundaries.
	 *********************************************************************/
	class element
	{
	public:
		virtual ~element () {}
		/*!*******************************************************************
		 * \brief This virtual function should update the element by one step
		 *********************************************************************/
		virtual void update () = 0;
	};
	
	class element_1D : public element
	{
	public:
		element_1D (int i_n) {
			int i;
			n = i_n;
			n_scalars = 0;
			
			cell.resize (i_n);
			for (i = 0; i < i_n; ++i) {
				cell [i] = i;
			}
		}
		
		virtual ~element_1D () {}
		
		inline void add_scalar (int index) {
			++n_scalars;
			
			scalars [index].resize (n);
		}
		
		double boundary_top (int deriv);
		
		double boundary_bottom (int deriv);
	
	protected:
		int n; //!< The number of elements in each 1D array
		int n_scalars;
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
	class diffusion_element : public element_1D
	{
	public:
		/*!*******************************************************************
		 * \param i_n The number of elements in each 1D data array
		 * \param i_flags Flags for the boundary conditions and evaluation
		 *********************************************************************/
		diffusion_element (int i_n, int i_flags);
		virtual ~diffusion_element () {}
		/*!*******************************************************************
		 * \brief This calculates diffusion and output for a constant timestep
		 *********************************************************************/
		void update ();
	
	private:
		int flags; //!< Flags for the boundary conditions and evaluation
				
		std::unique_ptr<diffusion::collocation_chebyshev_1D> diffusion_plan; //!< The diffusion implementation
		std::unique_ptr<io::incremental_output> angle_stream; //!< An implementation to output in angle space
		std::unique_ptr<io::simple_output> failsafe_dump; //!< An implementation to dump in case of failure
		fftw_plan fourier_plan; //!< The fft implementation
	};
} /* element */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
