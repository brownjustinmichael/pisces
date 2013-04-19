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

#include <memory>
#include "config.hpp"

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
	plan (int *i_flags_ptr = NULL) {
		default_flags = 0x00;
		if (!i_flags_ptr) {
			flags_ptr = &default_flags;
		} else {
			flags_ptr = i_flags_ptr;
		}
	}
	virtual ~plan () {}
	/*!*******************************************************************
	 * \brief Operate the plan on the data arrays contained in the class
	 * 
	 * The plan class serves as a wrapper for this function.
	 *********************************************************************/
	virtual void execute () = 0;
	
protected:
	int default_flags; //!< An integer set of default flags to use in case the user does not specify any flags
	int *flags_ptr;
};

class implicit_plan : public plan
{
public:
	implicit_plan (int i_n, double *i_matrix, int *i_flags_ptr = NULL) : plan (i_flags_ptr) {
		n = i_n;
		matrix = i_matrix;
	}
	virtual ~implicit_plan () {}

protected:
	int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
	double *matrix; //!< A double pointer to the input data
};

class explicit_plan : public plan
{
public:
	explicit_plan (int i_n, double *i_data_in, double *i_data_out, int *i_flags_ptr = NULL) : plan (i_flags_ptr) {
		n = i_n;
		data_in = i_data_in;
		data_out = i_data_out;
	}
	virtual ~explicit_plan () {}

protected:
	int n; //!< An integer number of data elements (grid points) that collocation_1D will be built to handle
	double *data_in; //!< A double pointer to the input data
	double *data_out; //!< A double pointer to the output data; if data_in == data_out, the operation is done in place (but inefficient)
};


#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
