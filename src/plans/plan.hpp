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
#include "io/parameters.hpp"
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
		datatype coeff;
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
		plan <datatype> (int *i_element_flags = NULL, int *i_component_flags = NULL, datatype i_coeff = 1.0) :
		coeff (i_coeff), 
		element_flags (i_element_flags),
		component_flags (i_component_flags) {}
		
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

			virtual const int type () const {
				return 0x00;
			}
			
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
			virtual std::shared_ptr <plan <datatype>> instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
					std::shared_ptr <plan <datatype>> result = _instance (grids, matrices, i_data_in, i_data_out, i_element_flags, i_component_flags);
					if (result) result->coeff *= coeff;
					return result;
				}

		protected:
			virtual std::shared_ptr <plan <datatype>> _instance (grids::grid <datatype> **grids, datatype **matrices, datatype *i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const = 0;
		};

		class factory_container
		{
		private:

		public:
			std::vector <std::shared_ptr <factory>> facts;

			factory_container (std::shared_ptr <factory> i_fact = NULL) {
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

			factory_container operator* (YAML::Node &node) {
				if (node.IsDefined ()) {
					return *this * node.as <datatype> ();
				}
				return factory_container ();
			}

			
		};
	};

	template <class datatype>
	typename plan <datatype>::factory_container operator* (typename plan <datatype>::factory_container i_container, datatype scalar) {
		typename plan <datatype>::factory_container container (i_container);
		for (int i = 0; i < (int) i_container.facts.size (); ++i)
		{
			container.facts [i]->coeff *= scalar;
		}
		return container;
	}

	template <class datatype>
	typename plan <datatype>::factory_container operator* (datatype scalar, typename plan <datatype>::factory_container i_container) {
		return i_container * scalar;
	}

	template <class datatype>
	datatype operator* (YAML::Node &node, datatype other) {
		return other * node;
	}

} /* plans */
#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
