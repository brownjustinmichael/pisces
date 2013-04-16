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
#include "../plan.hpp"
#include "../collocation/collocation.hpp"
#include "../io/io.hpp"
#include "../solver/solver.hpp"

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
		
		inline void add_scalar (int name) {
			++n_scalars;
			
			scalars [name].resize (n);
		}
		
		inline double index (int name, int i) {
			return scalars [name] [i];
		}
		
	protected:
		int n; //!< The number of elements in each 1D array
		int n_scalars;
		std::vector<int> cell; //!< An integer array for tracking each cell number for output
		std::map<int, std::vector<double>> scalars; //!< A vector of scalar vectors
		std::shared_ptr<collocation::chebyshev_grid> grid;
	};

	/*!*******************************************************************
	 * \brief A simple implementation of the element class with diffusion
	 * 
	 * This class contains a full element's capacity to run a single 
	 * element diffusion in 1D. It cannot yet tie to other elements or 
	 * calculate its own timestep. 
	 *********************************************************************/
	class diffusion_element_1D : public element_1D
	{
	public:
		/*!*******************************************************************
		 * \param i_n The number of elements in each 1D data array
		 * \param i_flags Flags for the boundary conditions and evaluation
		 *********************************************************************/
		diffusion_element_1D (int i_n, int i_flags);
		virtual ~diffusion_element_1D () {}
		/*!*******************************************************************
		 * \brief This calculates diffusion and output for a constant timestep
		 *********************************************************************/
		void update ();
		
		inline double derivative (int name, int deriv, int k) {
			if (deriv == 0) {
				return index (name, k);
			} else {
				int i;
				double derivative = 0.0;
				for (i = 0; i < n; ++i) {
					derivative += scalars [name] [i] * grid->index (deriv, i, k);
				}
				return derivative;
			}
		};
		
	private:
		int flags; //!< Flags for the boundary conditions and evaluation
		double timestep;
		double previous_timestep;
		std::vector<double> matrix;
		std::unique_ptr<plan> implicit_diffusion; //!< The diffusion implementation
		std::unique_ptr<plan> explicit_diffusion; //!< The diffusion implementation
		std::unique_ptr<solver::solver> matrix_solver;
		std::unique_ptr<io::output> angle_stream; //!< An implementation to output in angle space
		std::unique_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		std::unique_ptr<plan> fourier_transform; //!< The fft implementation
	};
} /* element */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
