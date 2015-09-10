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
#include "versions/version.hpp"
#include "logger/logger.hpp"
#include "io/parameters.hpp"
#include "grids/grid.hpp"
#include "grids/variable.hpp"
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


enum states
{
	real_real = 0,
	real_spectral = 1,
	spectral_spectral = 2
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
		datatype coeff;
		datatype *data_in;
		datatype *data_out;
		grids::variable <datatype> &var_in; //!< A datatype pointer to the input data
		grids::variable <datatype> &var_out; //!< A datatype pointer to the input data
		int &element_flags; //!< A pointer to the integer global flags
		int &component_flags; //!< A pointer to the integer local flags

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
		plan (grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, int state_in = 0, int state_out = 0, datatype i_coeff = 1.0) :
		coeff (i_coeff), 
		data_in (i_data_in.ptr (state_in)),
		data_out (i_data_out.ptr (state_out)),
		var_in (i_data_in),
		var_out (i_data_out),
		element_flags (i_data_in.element_flags),
		component_flags (i_data_in.component_flags) {}

		plan (grids::variable <datatype> &i_data_in, int state = 0, datatype i_coeff = 1.0) :
		plan (i_data_in, i_data_in, state, state, i_coeff) {}
		
		virtual ~plan () {}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("2.1.1.0");
			return version;
		}
		
		virtual void setup () = 0;
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () = 0;

		virtual bool implicit () {
			return false;
		}

		class factory_container;
		
		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce a plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory
		{
		public:
			datatype coeff;

			enum type {
				impl = 0x01,
				expl = 0x02,
				real = 0x04,
				solv = 0x08
			};

			factory (datatype i_coeff = 0.0) : coeff (i_coeff) {}

			virtual ~factory () {}

			/*!**********************************************************************
			 * \brief The abstract instance creating class
			 * 
			 * \param grids An array of grid objects that define the data grid
			 * \param matrices An array of the solver datatype matrices
			 * \param i_data_in The datatype array of the input data
			 * \param i_data_out The datatype array of the output data
			 * \param i_element_flags A pointer to the integer flags associated with the element on the whole
			 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
			 * 
			 * This method creates a shared_ptr to an implicit plan instance. The benefit to this inclusion is that the instance method can be called in a uniform way and hide communication of grid and matrix information from the user. If a plan would be created that would not do anything (e.g. something with a coefficient of 0.0), this will return a NULL shared pointer.
			 ************************************************************************/
			virtual std::shared_ptr <plan <datatype>> instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
					return _instance (matrices, i_data_in, i_data_out);
				}

		protected:
			virtual std::shared_ptr <plan <datatype>> _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const = 0;
		};

		class factory_container
		{
		private:

		public:
			std::vector <std::shared_ptr <factory>> facts;

			factory_container () {}

			factory_container (std::shared_ptr <factory> i_fact) {
				facts.push_back (i_fact);
			}

			virtual ~factory_container () {}

			factory_container operator+ (factory_container j_container) {
				factory_container container (*this);
				for (int i = 0; i < (int) j_container.facts.size (); ++i)
				{
					container.facts.push_back (j_container.facts [i]);
				}
				return container;
			}

			factory_container operator+ (std::shared_ptr <plan <datatype>::factory> i_factory) {
				return *this + factory_container (i_factory);
			}

			factory_container operator- (factory_container j_container) {
				return *this + (j_container * -1.0);
			}

			factory_container operator- (std::shared_ptr <plan <datatype>::factory> i_factory) {
				return *this - factory_container (i_factory);
			}

			factory_container operator+ (grids::variable <datatype> &var) {
				return *this + src (var);
			}

			factory_container operator+ (datatype scalar) {
				return *this + constant (scalar);
			}

			factory_container operator* (YAML::Node &node) {
				// if (node.IsDefined ()) {
					return *this * node.as <datatype> ();
				// }
				WARN ("Missing parameter... Assuming to be 0.0")
				return factory_container ();
			}

			
		};
	};

	plan <double>::factory_container operator+ (std::shared_ptr <plan <double>::factory> i_factory, std::shared_ptr <plan <double>::factory> j_factory);

	plan <double>::factory_container operator+ (std::shared_ptr <plan <double>::factory> i_factory, plan <double>::factory_container j_container);

	template <class type>
	typename plan <double>::factory_container operator- (std::shared_ptr <typename plan <double>::factory> i_factory, type i_other) {
		return i_factory + (-1.) * i_other;
	}

	template <class datatype>
	typename plan <datatype>::factory_container operator* (typename plan <datatype>::factory_container i_container, datatype scalar);

	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> operator* (std::shared_ptr <typename plan <datatype>::factory> i_factory, datatype scalar) {
		i_factory->coeff *= scalar;
		return i_factory;
	}

	plan <double>::factory_container operator- (std::shared_ptr <plan <double>::factory> i_factory);

	std::shared_ptr <plan <double>::factory> operator* (std::shared_ptr <plan <double>::factory> i_factory, YAML::Node node);

	template <class datatype>
	typename plan <datatype>::factory_container operator* (datatype scalar, typename plan <datatype>::factory_container i_container);

	template <class datatype>
	datatype operator* (YAML::Node node, datatype other) {
		return other * node;
	}

} /* plans */

template <class datatype>
std::shared_ptr <plans::plan <datatype>> NULL_plan () {
	return std::shared_ptr <plans::plan <datatype>> ();
}


#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
