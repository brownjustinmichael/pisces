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

#include "collocation.hpp"
#include "../config.hpp"

/*!*******************************************************************
 * \brief A set of indices to be used with the element scalars for convenience
 * 
 * Indices less than zero are "resettable." The null index is 
 * equivalent to a NULL pointer.
 *********************************************************************/
enum index {
	null = 00,
	
	position = 01, x_position = 01, x_pos = 01,
	y_position = 02, y_pos = 02, 
	z_position = 03, z_pos = 03,
	velocity = 11, vel = 11, x_velocity = 11, x_vel = 11,
	y_velocity = 12, y_vel = 12,
	z_velocity = 13, z_vel = 13,
	pressure = 20, pres = 20,
	temperature = 21, temp = 21,
	composition = 22, comp = 22,
	
	rhs = -01,
	vel_rhs = -11,
	temp_rhs = -21
};

/*!*******************************************************************
 * \brief A set of flags to be used with the plan class
 *********************************************************************/
enum plan_flags {
	unchanged_timestep = 0x400,
	transformed = 0x10
};

namespace bases
{
	template <class datatype>
	class element;
	
	class messenger;
	
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
		/*!**********************************************************************
		 * \param i_element_ptr A pointer to the associated element
		 ************************************************************************/
		plan (element <datatype>* i_element_ptr, int flags = 0x00);
		
		virtual ~plan () {}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () {
			if (!element_ptr) {
				ERROR ("Plan not initialized completely.");
				throw 0;
			}
		}
			
	protected:
		element <datatype>* element_ptr; //!< A pointer to the element with which the plan is associated
		int flags; // The integer plan execution flags
		int default_flags; //!< An integer set of default flags to use in case the user does not specify any flags
		int *flags_ptr; //!< A pointer to the integer element execution flags
		messenger* messenger_ptr; //!< A pointer to the messenger associated with the element
	};

	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output.
	 *********************************************************************/
	template <class datatype>
	class explicit_plan : public plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data_in The integer scalar index of the input
		 * \param i_data_out The integer scalar index of the output
		 * \copydoc plan::plan ()
		 *********************************************************************/
		explicit_plan (element <datatype>* i_element_ptr, int i_n, datatype* i_data_in, datatype* i_data_out = NULL, int flags = 0x00);
	
		virtual ~explicit_plan () {}
	
		/*!*******************************************************************
		 * \copydoc bases::plan::execute ()
		 *********************************************************************/
		virtual void execute ();
		
	protected:
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		datatype* data_in; //!< A datatype pointer to the input data
		datatype* data_out; //!< A datatype pointer to the output data
	};
	
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to implicit methods
	 * 
	 * These plans produce output in a square matrix but take no input.
	 *********************************************************************/
	template <class datatype>
	class implicit_plan : public plan <datatype>
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in a row of the square i_matrix
		 * \param i_grid A shared pointer to the collocation grid object
		 * \param i_matrix The datatype matrix to be updated
		 * \copydoc plan::plan ()
		 *********************************************************************/
		implicit_plan (element <datatype>* i_element_ptr, int i_n, collocation_grid <datatype>* i_grid, datatype *i_matrix, int i_flags = 0x00) : 
		plan <datatype> (i_element_ptr, i_flags), 
		n (i_n),
		grid (i_grid),
		matrix (i_matrix) {}

		virtual ~implicit_plan () {}
		
		/*!*******************************************************************
		 * \copydoc plan::execute ()
		 *********************************************************************/
		virtual void execute ();
		
	protected:
		int n; //!< An integer number of data elements
		collocation_grid <datatype>* grid; //!< A shared pointer to the grid
		datatype *matrix; //!< A datatype pointer to the input data
	};
} /* bases */

#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
