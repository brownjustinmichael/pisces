/*!***********************************************************************
 * \file bases/element.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef ELEMENT_HPP_IUTSU4TQ
#define ELEMENT_HPP_IUTSU4TQ

#include "plan.hpp"
#include "solver.hpp"
#include "../io/io.hpp"
#include "collocation.hpp"

/*!*******************************************************************
 * \brief A set of indices to be used with the element scalars for convenience
 *********************************************************************/
enum index {
	position = 00, x_position = 00, x_pos = 00,
	y_position = 01, y_pos = 01, 
	z_position = 02, z_pos = 02,
	velocity = 10, x_velocity = 10, x_vel = 10,
	y_velocity = 11, y_vel = 11,
	z_velocity = 12, z_vel = 12,
	pressure = 20, pres = 20,
	temperature = 21, temp = 21,
	composition = 22, comp = 22,
	rhs = 30
};

namespace bases
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
		/*!*******************************************************************
		* \param i_flags An integer set of execution flags
		*********************************************************************/
		element (int i_flags) {
			flags = i_flags;
		
			n_explicit_grid_plans = 0;
			n_explicit_space_plans = 0;
			n_implicit_plans = 0;
			n_boundaries = 0;
		
			previous_timestep = 0.0;
		}
		
		virtual ~element () {}
		/*!*******************************************************************
		 * \brief Calculate the matrix terms to be used in update
		 *********************************************************************/
		virtual void calculate ();
	
		/*!*******************************************************************
		 * \brief Execute the boundary conditions
		 *********************************************************************/
		virtual void execute_boundaries ();
	
		/*!*******************************************************************
		 * \brief Update the element
		 *********************************************************************/
		virtual void update ();
		
		/*!*******************************************************************
		 * \brief Add a scalar to the element
		 * 
		 * \param name The integer name from the index enumeration
		 *********************************************************************/
		virtual void add_scalar (int name) = 0;
		
		/*!*******************************************************************
		 * \brief Get the double pointer to the named scalar
		 * 
		 * \param name The integer name from the index enumeration
		 * 
		 * \return A double pointer to the first element of the named scalar
		 *********************************************************************/
		virtual double &operator[] (int name) = 0;
		
		virtual double &operator () (int name, int index = 0) {
			return (&((*this) [name])) [index];
		}
	
		/*!*******************************************************************
		 * \brief Adds an explicit plan to be evaluated in grid space
		 * 
		 * \param i_plan A unique pointer to the plan to add
		 *********************************************************************/
		inline void add_explicit_grid_plan (std::unique_ptr<plan> i_plan) {
			++n_explicit_grid_plans;
			explicit_grid_plans.push_back (std::move (i_plan));
		}

		/*!*******************************************************************
		 * \brief Adds an explicit plan to be evaluated in normal space
		 * 
		 * \param i_plan A unique pointer to the plan to add
		 *********************************************************************/	
		inline void add_explicit_space_plan (std::unique_ptr<plan> i_plan) {
			++n_explicit_space_plans;
			explicit_space_plans.push_back (std::move (i_plan));
		}
	
		/*!*******************************************************************
		 * \brief Adds an implicit plan to be evaluated
		 * 
		 * \param i_plan A unique pointer to the plan to add
		 *********************************************************************/
		inline void add_implicit_plan (std::unique_ptr<plan> i_plan) {
			++n_implicit_plans;
			implicit_plans.push_back (std::move (i_plan));
		}
	
		/*!*******************************************************************
		 * \brief Adds a boundary condition to execute in normal space
		 * 
		 * \param i_plan A unique pointer to the plan to add
		 *********************************************************************/
		inline void add_boundary (std::unique_ptr<plan> i_plan) {
			++n_boundaries;
			boundaries.push_back (std::move (i_plan));
		}
	
protected:
	int flags; //!< An integer set of execution flags
	int n_explicit_grid_plans; //!< The number of explicit grid plans to execute
	int n_explicit_space_plans; //!< The number of explicit space plans to execute
	int n_implicit_plans; //!< The number of implicit plans to execute
	int n_boundaries; //!< The number of boundary conditions to execute
	
	double timestep; //!< The double timestep length
	double previous_timestep; //!< The double duration of the previous timestep

	std::shared_ptr<collocation_grid> grid; //!< A shared pointer to the collocation grid
	
	std::unique_ptr<plan> transform_forward; //!< A unique pointer to the forward transform
	std::unique_ptr<solver> matrix_solver; //!< A unique pointer to the matrix solver
	std::unique_ptr<io::output> angle_stream; //!< An implementation to output in angle space
	std::unique_ptr<io::output> failsafe_dump; //!< An implementation to dump in case of failure
		
	std::vector<std::unique_ptr<plan>> explicit_grid_plans; //!< A vector of unique pointers to explicit grid plans to be executed
	std::vector<std::unique_ptr<plan>> explicit_space_plans; //!< A vector of unique pointers to explicit space plans to be executed
	std::vector<std::unique_ptr<plan>> implicit_plans; //!< A vector of unique pointers to implicit plans to be executed
	std::vector<std::unique_ptr<plan>> boundaries; //!< A vector of unique pointers to boundary conditions to be executed
};
} /* bases */

#endif /* end of include guard: ELEMENT_HPP_IUTSU4TQ */
