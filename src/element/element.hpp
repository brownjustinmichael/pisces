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
#include <algorithm>
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
		virtual void calculate () = 0;
		virtual void execute_boundaries () = 0;
		virtual void update () = 0;
	};
	
	class element_1D : public element
	{
	public:
		element_1D (int i_n, int i_flags) {
			int i;
			n = i_n;
			flags = i_flags;
			n_explicit_grid_plans = 0;
			n_explicit_space_plans = 0;
			n_implicit_plans = 0;
			
			cell.resize (i_n);
			for (i = 0; i < i_n; ++i) {
				cell [i] = i;
			}
		}
		
		virtual ~element_1D () {}
		
		inline void add_scalar (int name) {
			scalars [name].resize (n + 1);
		}
		
		inline double index (int name, int i) {
			return scalars [name] [i];
		}
		
		inline void add_explicit_grid_plan (std::unique_ptr<plan> i_plan) {
			++n_explicit_grid_plans;
			explicit_grid_plans.push_back (std::move (i_plan));
		}
		
		inline void add_explicit_space_plan (std::unique_ptr<plan> i_plan) {
			++n_explicit_space_plans;
			explicit_space_plans.push_back (std::move (i_plan));
		}
		
		inline void add_implicit_plan (std::unique_ptr<plan> i_plan) {
			++n_implicit_plans;
			implicit_plans.push_back (std::move (i_plan));
		}
		
		inline void add_boundary (std::unique_ptr<plan> i_plan) {
			++n_boundaries;
			boundaries.push_back (std::move (i_plan));
		}
		
		virtual void calculate ();
		virtual void execute_boundaries ();
		virtual void update ();
		
	protected:
		int n; //!< The number of elements in each 1D array
		int flags;
		int n_explicit_grid_plans;
		int n_explicit_space_plans;
		int n_implicit_plans;
		int n_boundaries;
		double timestep;
		double previous_timestep;
		std::vector<int> cell; //!< An integer array for tracking each cell number for output
		std::map<int, std::vector<double>> scalars; //!< A vector of scalar vectors
		std::vector<std::unique_ptr<plan>> explicit_grid_plans;
		std::vector<std::unique_ptr<plan>> explicit_space_plans;
		std::vector<std::unique_ptr<plan>> implicit_plans;
		std::vector<std::unique_ptr<plan>> boundaries;
		std::unique_ptr<plan> transform_forward;
		std::unique_ptr<solver::solver> matrix_solver;
		std::shared_ptr<collocation::chebyshev_grid> grid;
		std::unique_ptr<io::output> angle_stream; //!< An implementation to output in angle space
		std::unique_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
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
		std::vector<double> matrix;
		std::unique_ptr<plan> implicit_diffusion; //!< The diffusion implementation
		std::unique_ptr<plan> explicit_diffusion; //!< The diffusion implementation
	};
} /* element */

#endif /* end of include guard: ELEMENT_HPP_3SURDTOH */
