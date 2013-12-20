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

/*!*******************************************************************
 * \brief A set of indices to be used with the element scalars for convenience
 * 
 * Indices less than zero are "resettable." The null index is 
 * equivalent to a NULL pointer.
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
	
	rhs = -01,
	x_vel_explicit_rhs = -11,
	x_vel_real_rhs = -12,
	x_vel_implicit_rhs = -13,
	vel_explicit_rhs = -16, z_vel_explicit_rhs = -16,
	vel_real_rhs = -17, z_vel_real_rhs = -17,
	vel_implicit_rhs = -18, z_vel_implicit_rhs = -18,
	temp_explicit_rhs = -21,
	temp_implicit_rhs = -22,
	temp_real_rhs = -23,
	stream_explicit_rhs = -31,
	stream_implicit_rhs = -32,
	stream_real_rhs = -33
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
	only_forward_horizontal = 0x08
};

enum transform_flags {
	forward_horizontal = 0x01,
	forward_vertical = 0x02,
	inverse_horizontal = 0x04,
	inverse_vertical = 0x08
}

namespace bases
{
	typedef std::map <int, int> flags;
	
	/*!*******************************************************************
	* \brief The basic functional unit, containing a recipe for execution
	* 
	* An implemented plan class contains the operator and the addresses of 
	* all the relevant data arrays to operate on. Each plan need be 
	* constructed only once and can run any number of times each timestep.
	*********************************************************************/
	template <class datatype>
	class plan
	{
	public:
		virtual ~plan () {
			// printf ("Destroying bases plan\n");
		}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute (bases::flags &element_flags) = 0;
	};
} /* bases */

template <class datatype>
class change_flag : public bases::plan <datatype>
{
public:
	change_flag (int i_flag, bool i_boolean) :
	flag (i_flag),
	boolean (i_boolean) {}
	
	virtual ~change_flag () {
		printf ("Destroying change flag\n");
	}
	
	virtual void execute (bases::flags &element_flags) {
		if (boolean) {
			element_flags [state] |= flag;
		} else {
			element_flags [state] &= ~flag;
		}
	}

private:
	int flag;
	bool boolean;
};
#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
