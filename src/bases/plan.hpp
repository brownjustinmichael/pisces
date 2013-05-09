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

#include <memory>
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
	rhsd2 = -02
};

enum plan_flags {
	implicit_started = 0x100,
	explicit_started = 0x200
};

namespace bases
{
	class element;
	
	class boundary;

	/*!*******************************************************************
	* \brief The basic functional unit, containing a recipe for execution
	* 
	* An implemented plan class contains the operator and the addresses of 
	* all the relevant data arrays to operate on. Each plan need be 
	* constructed only once and can run any number of times each timestep.
	*********************************************************************/
	class plan
	{
	public:
		/*!*******************************************************************
		 * \param i_flags_ptr A pointer to the integer execution flags
		 *********************************************************************/
		plan () {
			MTRACE ("Instantiating...");
			element_ptr = NULL;
			MTRACE ("Instantiated.");
		}
		
		virtual ~plan () {}
		
		virtual void associate (element* i_element_ptr);
		
		virtual void setup_boundary (boundary* i_boundary) {}
		
		virtual void boundary (boundary* i_boundary) {}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () {
			if (!element_ptr) {
				MERROR ("Plan not initialized completely.");
				throw 0;
			}
		}
			
	protected:
		int logger;
		int default_flags; //!< An integer set of default flags to use in case the user does not specify any flags
		int *flags_ptr; //!< A pointer to the integer execution flags
		element* element_ptr;
	};

	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output.
	 *********************************************************************/
	class explicit_plan : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in the data
		 * \param i_data_in The double vector of input
		 * \param i_data_out The double vector of output (if NULL, use i_data_in)
		 * \copydoc plan::plan ()
		 *********************************************************************/
		explicit_plan (int i_n, int i_name_in, int i_name_out = null) : plan () {
			MTRACE ("Instantiating...");
			n = i_n;
			name_in = i_name_in;
			if (!i_name_out) {
				name_out = i_name_in;
			} else {
				name_out = i_name_out;
			}
			MTRACE ("Instantiated.");
		}
	
		virtual ~explicit_plan () {}
	
		virtual void associate (element* i_element_ptr);
		
		virtual void execute ();
		
	protected:
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		int name_in; //!< A double pointer to the input data
		int name_out; //!< A double pointer to the output data
		double* data_in;
		double* data_out;
	};
	
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to implicit methods
	 * 
	 * These plans produce output in a square matrix but take no input.
	 *********************************************************************/
	class implicit_plan : public plan
	{
	public:
		/*!*******************************************************************
		 * \param i_n The integer number of elements in a row of the square i_matrix
		 * \param i_matrix The double matrix to be updated
		 * 
		 * \copydoc plan::plan ()
		 *********************************************************************/
		implicit_plan (int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix) : plan () {
			n = i_n;
			grid = i_grid;
			matrix = i_matrix;
		}

		virtual ~implicit_plan () {}
		
		virtual void execute ();
		
	protected:
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		std::shared_ptr <bases::collocation_grid> grid;
		double *matrix; //!< A double pointer to the input data
	};
} /* bases */

#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
