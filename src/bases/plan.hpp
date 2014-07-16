/*!***********************************************************************
 * \file bases/plan.hpp
 * Spectral Element
 * 
 * This file provides the abstract base class plan, from which all code 
 * executions should derive.
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PLAN_HPP_S9YPWHOM
#define PLAN_HPP_S9YPWHOM

#include "../config.hpp"
#include <map>

#define NANTRACK

/*!*******************************************************************
 * \brief A set of indices to be used with the element scalars for convenience
 * 
 * Indices less than zero are "resettable." The null index is equivalent to a NULL pointer.
 *********************************************************************/
enum index {
	null = 0x0,
	state = 0x0,
	
	x_position = 01, x_pos = 01,
	y_position = 02, y_pos = 02, 
	position = 03, z_position = 03, z_pos = 03,
	vel = 11, x_velocity = 11, x_vel = 11,
	y_velocity = 12, y_vel = 12,
	velocity = 13, z_velocity = 13, z_vel = 13,
	stream = 14, streamfunction = 14,
	pressure = 20, pres = 20,
	temperature = 21, temp = 21,
	composition = 22, comp = 22,
	test = 100, test_x = 101, text_z = 103
};

/*!*******************************************************************
 * \brief A set of flags to be used with the plan class
 *********************************************************************/
enum plan_flags {
	unchanged_timestep = 0x400,
	implicit_set = 0x4000,
	transformed_horizontal = 0x10,
	transformed_vertical = 0x8000,
	no_transform = 0x04,
	only_forward_horizontal = 0x08,
	changed = 0x10000,
	transform_output = 0x20000,
	normal_output = 0x40000,
	profile_only = 0x80000,
	timestep_only = 0x100000,
	no_solve = 0x200000,
	solved = 0x400000
};

enum solve_element_flags {
	x_solve = 0x20,
	z_solve = 0x80
};

namespace bases
{
	/*!*******************************************************************
	* \brief The basic functional unit, containing a recipe for execution
	* 
	* An implemented plan class contains the operator and the addresses of all the relevant data arrays to operate on. Each plan need be constructed only once and can run any number of times each timestep.
	*********************************************************************/
	template <class datatype>
	class plan
	{
	protected:
		int *element_flags; //!< A pointer to the integer global flags
		int *component_flags; //!< A pointer to the integer local flags
		
	public:
		/*!**********************************************************************
		* \param i_element_flags A pointer to the integer global flags
		* \param i_component_flags A pointer to the integer local flags
		 ************************************************************************/
		plan <datatype> (int *i_element_flags = NULL, int *i_component_flags = NULL) :
		element_flags (i_element_flags),
		component_flags (i_component_flags)  {}
		
		virtual ~plan () {}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () = 0;
	};
} /* bases */
#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
