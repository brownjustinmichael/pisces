/*!***********************************************************************
 * \file plan.hpp
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

#include "versions/version.hpp"
#include "logger/logger.hpp"
#include "grids/grid.hpp"
#include <map>

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
	solved = 0x4000000,
	plans_setup = 0x8000000
};

/*!**********************************************************************
 * \brief A set of flags to determine the direction of the solve
 ************************************************************************/
enum solve_element_flags {
	x_solve = 0x20,
	z_solve = 0x80
};

/*!**********************************************************************
 * \namespace plans
 * 
 * \brief A namespace containing the plans that do the actual operations on data.
 ************************************************************************/
namespace plans
{
	template <class datatype>
	struct variable	{
		std::string name;
		int *element_flags;
		int *component_flags;
		std::vector <grids::grid <datatype> *> grid_ptrs;	
	};	

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
		enum types {
			pre = 0x01,
			mid = 0x02,
			post = 0x04,
			pre_solve = 0x08
		};

		/*!**********************************************************************
		* \param i_element_flags A pointer to the integer global flags
		* \param i_component_flags A pointer to the integer local flags
		 ************************************************************************/
		plan <datatype> (int *i_element_flags = NULL, int *i_component_flags = NULL) :
		element_flags (i_element_flags),
		component_flags (i_component_flags)  {}
		
		virtual ~plan () {}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("1.1.1.0");
			return version;
		}
		
		virtual void setup () = 0;
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () = 0;

		virtual int type () = 0;
	};
} /* plans */
#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
