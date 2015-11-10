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

/**
 * @brief A list of states to be used as the states for the variable class
 */
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
	/*!*******************************************************************
	* \brief The basic functional unit, containing a recipe for execution
	* 
	* An implemented plan class contains the operator and the addresses of all the relevant data arrays to operate on. Each plan need be constructed only once and can run any number of times each timestep.
	*********************************************************************/
	template <class datatype>
	class plan
	{
	protected:
		datatype coeff; //!< A coefficient to multiply the results of the plan operation
		datatype *data_in; //!< A pointer to the data to operate on
		datatype *data_out; //!< A pointer to the location to output to
		grids::variable &var_in; //!< A datatype pointer to the input data
		grids::variable &var_out; //!< A datatype pointer to the input data
		int &element_flags; //!< A pointer to the integer global flags
		int &component_flags; //!< A pointer to the integer local flags

	public:
		/**
		 * @brief A set of types of plans to suggest to the class which inputs and outputs to use
		 */
		enum types {
			pre = 0x01,
			mid = 0x02,
			post = 0x04,
			pre_solve = 0x08
		};

		/*!**********************************************************************
		 * \param i_data_in A reference to the input data
		 * @param i_data_out A reference to the output data
		 * @param state_in An integer of which state of data_in to use as the input
		 * @param state_out An integer of which state of data_out to use as the output
		 * @param i_coeff A coefficient to multiply the results of the plan operation
		 ************************************************************************/
		plan (grids::variable &i_data_in, grids::variable &i_data_out, int state_in = 0, int state_out = 0, datatype i_coeff = 1.0) :
		coeff (i_coeff), 
		data_in (i_data_in.ptr (state_in)),
		data_out (i_data_out.ptr (state_out)),
		var_in (i_data_in),
		var_out (i_data_out),
		element_flags (i_data_in.element_flags),
		component_flags (i_data_in.component_flags) {}

		/*!**********************************************************************
		 * @brief An in-place plan operation where the output will overwrite the input
		 * 
		 * \param i_data_in A reference to the input data
		 * @param state An integer of which state of data_in to use as the input and output
		 * @param i_coeff A coefficient to multiply the results of the plan operation
		 ************************************************************************/
		plan (grids::variable &i_data_in, int state = 0, datatype i_coeff = 1.0) :
		plan (i_data_in, i_data_in, state, state, i_coeff) {}
		
		virtual ~plan () {}
		
		/*!**********************************************************************
		 * \return The version of the class
		 ************************************************************************/
		static versions::version& version () {
			static versions::version version ("2.1.1.0");
			return version;
		}
		
		/**
		 * @brief Set up any auxiliary information needed by the plan
		 * @details Defaults to an empty function for a new plan, but can be used to set up the plan again in case the implicit matrices become invalid.
		 */
		virtual void setup () {}
		
		/*!*******************************************************************
		* \brief Operate the plan on the data arrays contained in the class
		* 
		* The plan class serves as a wrapper for this function.
		*********************************************************************/
		virtual void execute () = 0;

		/**
		 * @brief Returns whether the plan is implicit
		 * 
		 * @return True if the plan is implicit (i.e. has a matrix to factorize), false otherwise
		 */
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
			datatype coeff; //!< The coefficient to be applied to the generated plan

			/**
			 * @brief A factory type used to identify which kind of plan we could be dealing with
			 */
			enum type {
				impl = 0x01,
				expl = 0x02,
				real = 0x04,
				solv = 0x08
			};

			/**
			 * @param i_coeff The coefficient to be applied to the generated plan
			 */
			factory (datatype i_coeff = 0.0) : coeff (i_coeff) {}

			virtual ~factory () {}

			/*!**********************************************************************
			 * \brief Create an instance of the plan from the factory
			 * 
			 * \param matrices An array of the solver datatype matrices
			 * \param i_data_in A reference to the variable input data
			 * \param i_data_out A reference to the variable output data
			 * 
			 * This method creates a shared_ptr to an implicit plan instance. The benefit to this inclusion is that the instance method can be called in a uniform way and hide communication of grid and matrix information from the user. If a plan would be created that would not do anything (e.g. something with a coefficient of 0.0), this will return a NULL shared pointer.
			 ************************************************************************/
			virtual std::shared_ptr <plan <datatype>> instance (datatype **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const {
					return _instance (matrices, i_data_in, i_data_out);
				}

		protected:
			/*!**********************************************************************
			 * \brief The abstract instance creating method
			 * 
			 * \param matrices An array of the solver datatype matrices
			 * \param i_data_in A reference to the variable input data
			 * \param i_data_out A reference to the variable output data
			 * 
			 * This method should be overloaded for each plan to allow for convenient plan generation
			 ************************************************************************/
			virtual std::shared_ptr <plan <datatype>> _instance (datatype **matrices, grids::variable &i_data_in, grids::variable &i_data_out) const = 0;
		};

		/**
		 * @brief A container class that functions to hold many factories
		 * @details This class is designed to contain a series of factories in order to make plan instantiation easier and more object oriented.
		 */
		class factory_container
		{
		private:

		public:
			std::vector <std::shared_ptr <factory>> facts; //!< A vector of shared pointers to factory instances

			factory_container () {}

			/**
			 * @param i_fact A shared pointer to a factory as the first element of the container
			 */
			factory_container (std::shared_ptr <factory> i_fact) {
				facts.push_back (i_fact);
			}

			virtual ~factory_container () {}

			/**
			 * @brief Returns a new container with the contents of both this one and another one
			 * @details This returns a new container with the contents of both this one and another one
			 * 
			 * @param j_container The container to add this one to
			 * @return A new container containing the contents of both containers
			 */
			factory_container operator+ (factory_container j_container) {
				factory_container container (*this);
				for (int i = 0; i < (int) j_container.facts.size (); ++i)
				{
					container.facts.push_back (j_container.facts [i]);
				}
				return container;
			}

			/**
			 * @brief Add a new factory to this container
			 * 
			 * @param i_factory The new factory to append to the container
			 * @return A new factory container containing the old contents and the new factory
			 */
			factory_container operator+ (std::shared_ptr <plan <datatype>::factory> i_factory) {
				return *this + factory_container (i_factory);
			}

			/**
			 * @brief Add a source term to this container
			 * @details If the added object is a variable, assume that a source term with that variable as the input should be added.
			 * 
			 * @param var A reference to the variable source
			 * @return A new factory continer with the old contents and the new source factory
			 */
			factory_container operator+ (grids::variable &var);

			/**
			 * @brief Add a constant term to this container
			 * @details If the added object is a simple scalar, add a constant source term
			 * 
			 * @param scalar The constant scalar to use as the source term
			 * @return A new factory container with the old contents and the new constant source factory
			 */
			factory_container operator+ (datatype scalar);

			/**
			 * @brief Append a new factory container with -1 times the coefficients contained
			 * @details Factory containers are designed to give the illusion of equations, so a subtraction actually represents the inclusion of another plan of opposite sign.
			 * 
			 * @param j_container The factory container to append to this one
			 * @return A new factory container with the old contents of this one and -1 times the contents of j_container
			 */
			factory_container operator- (factory_container j_container) {
				return *this + (j_container * -1.0);
			}

			/**
			 * @brief Append a new factory instance with -1 times the factory coefficient
			 * @details Factory containers are designed to give the illusion of equations, so a subtraction actually represents the inclusion of another plan of opposite sign.
			 * 
			 * @param i_factory The factory to append to this one
			 * @return A new factory container with the old contents of this one and -1 times the new factory
			 */
			factory_container operator- (std::shared_ptr <plan <datatype>::factory> i_factory) {
				return *this - factory_container (i_factory);
			}

			/**
			 * @brief Scale the contents of the factory container by a scalar
			 * @details Scale the coefficients of each factory contained by this container by the value of the scalar.
			 * 
			 * @param scalar The value by which to scale the contained coefficients
			 * @return A new factory container with the old contents scaled by scalar
			 */
			factory_container operator* (datatype scalar) {
				factory_container container (*this);
				for (int i = 0; i < (int) this->facts.size (); ++i)
				{
					container.facts [i]->coeff *= scalar;
				}
				return container;
			}

			/**
			 * @brief Scale the contents of the factory by a YAML::Node
			 * @details Assume that the node has the appropriate datatype (i.e. double or float) and scale the coefficients of the contained factories by that value. If the Node doesn't exist or has the wrong type, this will raise an error.
			 * 
			 * @param node A reference to the YAML Node by which to scale the contents
			 * @return A factory container with the old contents scaled by the node value
			 */
			factory_container operator* (YAML::Node &node) {
				// if (node.IsDefined ()) {
					return *this * node.as <datatype> ();
				// }
				WARN ("Missing parameter... Assuming to be 0.0")
				return factory_container ();
			}

			
		};
	};

	std::shared_ptr <typename plan <double>::factory> src (grids::variable &data_source, bool dealias);
	std::shared_ptr <typename plan <double>::factory> constant (double coeff);

	template <class datatype>
	typename plan <datatype>::factory_container plan <datatype>::factory_container::operator+ (grids::variable &var) {
		return *this + src (var, false);
	}

	template <class datatype>
	typename plan <datatype>::factory_container plan <datatype>::factory_container::operator+ (datatype scalar) {
		return *this + constant (scalar);
	}

	/**
	 * @brief Two shared pointers of factories add to a factory container containing both
	 * 
	 * @param i_factory A shared pointer to the first plan factory
	 * @param j_factory A shared pointer to the second plan factory
	 * @return The factory container containing both of the given plans
	 */
	plan <double>::factory_container operator+ (std::shared_ptr <plan <double>::factory> i_factory, std::shared_ptr <plan <double>::factory> j_factory);

	/**
	 * @brief A shared pointer of a factory and a factory container add to a factory container containing both
	 * 
	 * @param i_factory A shared pointer to a plan factory
	 * @param j_container A factory container
	 * @return The factory container containing the given plan factory and the contents of the container
	 */
	plan <double>::factory_container operator+ (std::shared_ptr <plan <double>::factory> i_factory, plan <double>::factory_container j_container);

	/**
	 * @brief Reworks subtraction of a shared pointer to a factory in terms of addition and multiplication
	 * 
	 * @param i_factory [description]
	 * @param i_other [description]
	 * 
	 * @return [description]
	 */
	template <class type>
	typename plan <double>::factory_container operator- (std::shared_ptr <typename plan <double>::factory> i_factory, type i_other) {
		return i_factory + (-1.) * i_other;
	}

	/**
	 * @brief Take the negative of the shared pointer of a factory (multiply the coefficient by -1)
	 * 
	 * @param i_factory The factory to be inverted
	 * @return A factory container containing the factory
	 */
	plan <double>::factory_container operator- (std::shared_ptr <plan <double>::factory> i_factory);

	/**
	 * @brief Multiply the coefficient of a shared pointer of a factory by a scalar
	 * 
	 * @param i_factory The shared pointer of the factory
	 * @param scalar The scalar to multiply the factory by
	 * 
	 * @return The shared pointer to the factory
	 */
	template <class datatype>
	std::shared_ptr <typename plan <datatype>::factory> operator* (std::shared_ptr <typename plan <datatype>::factory> i_factory, datatype scalar) {
		i_factory->coeff *= scalar;
		return i_factory;
	}

	/**
	 * @brief Multiply the coefficient of a shared pointer of a factory by a YAML Node
	 * @details This does first convert the YAML Node into a desired datatype, which currently must be specified in plan.cpp
	 * 
	 * @param i_factory The shared pointer of the factory
	 * @param node The YAML Node to multiply the coefficient of the factory by
	 * 
	 * @return The shared pointer to the factory
	 */
	std::shared_ptr <plan <double>::factory> operator* (std::shared_ptr <plan <double>::factory> i_factory, YAML::Node node);

	/**
	 * @brief If ever a plain datatype multiplies a factory container, reverse them
	 * 
	 * @param scalar The scalar value in the multiplication
	 * @param i_container The container object for which the coefficients should be multiplied
	 * 
	 * @return The factory container with the new coefficients
	 */
	template <class datatype>
	typename plan <datatype>::factory_container operator* (datatype scalar, typename plan <datatype>::factory_container i_container);

	/**
	 * @brief If ever in an instance where the YAML node is first in the multiplication, reverse them
	 * 
	 * @param node The YAML Node to reverse in the multiplication
	 * @param other Any other type
	 * 
	 * @return Something in the type of other
	 */
	template <class datatype>
	datatype operator* (YAML::Node node, datatype other) {
		return other * node;
	}

} /* plans */

/**
 * @return A null plan for convenience
 */
template <class datatype>
std::shared_ptr <plans::plan <datatype>> NULL_plan () {
	return std::shared_ptr <plans::plan <datatype>> ();
}

#endif /* end of include guard: PLAN_HPP_S9YPWHOM */
