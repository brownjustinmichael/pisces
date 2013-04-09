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
#include <fftw3.h>
#include "../diffusion/diffusion.hpp"
#include "../io/io.hpp"

namespace element
{
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

	/*!*******************************************************************
	 * \brief A simple implementation of the element class with diffusion
	 * 
	 * This class contains a full element's capacity to run a single 
	 * element diffusion in 1D. It cannot yet tie to other elements or 
	 * calculate its own timestep. 
	 *********************************************************************/
	class diffusion_element : public element
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
		int n; //!< The number of elements in each 1D array
		int flags; //!< Flags for the boundary conditions and evaluation
		
		std::vector<int> cell; //!< An integer array for tracking each cell number for output
		std::vector<double> position; //!< A double array containing position information
		std::vector<double> velocity; //!< A double array containing velocity information
		
		std::unique_ptr<diffusion::cheb_1D> diffusion_plan; //!< The diffusion implementation
		std::unique_ptr<io::incremental_output> angle_stream; //!< An implementation to output in angle space
		std::unique_ptr<io::simple_output> failsafe_dump; //!< An implementation to dump in case of failure
		fftw_plan fourier_plan; //!< The fft implementation
	};
} /* element */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
