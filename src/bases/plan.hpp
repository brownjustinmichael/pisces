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
 * \brief Function from BLAS that copies a double array to another in place
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 * \param dy The array to which the data are copied
 * \param incy A pointer to an integer spacing of elements in y
 *********************************************************************/
extern "C" void dcopy_ (int *n, double *dx, int *incx, double *dy, int *incy);

/*!*******************************************************************
 * \brief Function from BLAS that scales a double array
 * 
 * \param n A pointer to an integer number of elements in x to copy to y
 * \param da The double by which to scale the data
 * \param dx The array from which the data are copied
 * \param incx A pointer to an integer spacing of elements in x
 *********************************************************************/
extern "C" void dscal_ (int *n, double *da, double *dx, int *incx);

namespace bases
{
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
		plan (int *i_flags_ptr = NULL, int i_logger = -1) {
			logger = i_logger;
			TRACE (logger, "Instantiating...");
			default_flags = 0x00;
			if (!i_flags_ptr) {
				flags_ptr = &default_flags;
			} else {
				flags_ptr = i_flags_ptr;
			}
			TRACE (logger, "Instantiated.");
		}
		virtual ~plan () {}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () {}
			
	protected:
		int logger;
		int default_flags; //!< An integer set of default flags to use in case the user does not specify any flags
		int *flags_ptr; //!< A pointer to the integer execution flags
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
		explicit_plan (int i_n, double *i_data_in, double *i_data_out = NULL, int *i_flags_ptr = NULL, int i_logger = -1) : plan (i_flags_ptr, i_logger) {
			TRACE (logger, "Instantiating...");
			n = i_n;
			data_in = i_data_in;
			if (!i_data_out) {
				data_out = i_data_in;
			} else {
				data_out = i_data_out;
			}
			TRACE (logger, "Instantiated.");
		}
	
	virtual ~explicit_plan () {}
	
	virtual void execute () {
		plan::execute ();
	}

	protected:
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		double *data_in; //!< A double pointer to the input data
		double *data_out; //!< A double pointer to the output data
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
		implicit_plan (int i_n, std::shared_ptr<bases::collocation_grid> i_grid, double *i_matrix, int *i_flags_ptr = NULL, int i_logger = -1) : plan (i_flags_ptr, i_logger) {
			n = i_n;
			grid = i_grid;
			matrix = i_matrix;
		}

		virtual ~implicit_plan () {}
		
		virtual void execute () {
			plan::execute ();
		}
		

	protected:
		int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
		std::shared_ptr <bases::collocation_grid> grid;
		double *matrix; //!< A double pointer to the input data
	};
} /* bases */

#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
