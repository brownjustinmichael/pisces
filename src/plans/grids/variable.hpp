/*!***********************************************************************
 * \file variable.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VARIABLE_H__
#define VARIABLE_H__

#include "grid.hpp"

namespace grids
{
	/**
	 * @brief A set of flags for use with the variable class
	 */
	enum variable_flags
	{
		updated = 0x80000000
	};

	/**
	 * @brief A class designed to contain the values and state of a particular variable
	 * @details The variable class is allowed to have a number of different states (e.g. real-real, real-spectral), which means that this contains space for a few times the amount of information nominally present.
	 */
	class variable
	{
	protected:
		static std::vector <std::shared_ptr <variable>> tmps; //!< A place to store freely floating variables if necessary

		std::vector <grids::grid *> grids; //!< A vector of pointers to the grid objects associated with the system
		int dimensions; //!< Currently unimplemented, I'd like to allow variables to also be vectors on occasion
		std::vector <double> data; //!< The actual data contained
		int states; //!< The number of states in the variable
		// TODO This won't work. We need a way to keep references or entire variables
		std::vector <variable *> inner; //!< Some variables are compound, allowing them to be constructed from arithmetic operations on other variables, pointers to which are kept here
		std::vector <double *> vars; //!< Some variables are compound, allowing them to be constructed from arithmetic operations on other variables, pointers to the data of which are kept here
		std::vector <int> ops; //!< Some variables are compound, allowing them to be constructed from arithmetic operations on other variables, the arithmetic operations are kept here
		std::vector <int> inner_states; //!< A record of the previous states for each inner variable
		int ld; //!< The leading dimension of the system in the vertical direction
		int total; //!< The total number of data points in the variable (including states but not dimensions)

	public:
		int state; //!< The current state of the system (bad name; this is the current version)
		int last_update; //!< Which state was updated last
		int component_flags; //!< The integer flags associated with the variable
		int &element_flags; //!< A reference to the global element flags
		std::string name; //!< The name of the variable

		/**
		 * @brief Integer forms for the possible operators that can be used in compound variables
		 */
		enum operators {
			add = 0x01,
			sub = 0x02,
			mul = 0x03,
			div = 0x04
		};

		/**
		 * @brief Store a variable in the global temporary variables
		 * @details In case a variable is constructed that should be saved but that isn't associated with a particular name or data structure, it can be stored here.
		 * 
		 * @param other A reference to the variable to store
		 */
		static void store_var (std::shared_ptr <variable> &other) {
			tmps.push_back (other);
		}

		/**
		 * @brief If the stored variables need to be updated, this does so
		 * @details This calls update() on any stored variables.
		 */
		static void update_tmps () {
			for (auto iter = tmps.begin (); iter != tmps.end (); ++iter)
			{
				(*iter)->update ();
			}
		}

		/**
		 * @param i_grid_m A reference to the vertical grid
		 * @param i_element_flags A reference to the global element flags
		 * @param i_name The string name of the variable
		 * @param i_states The number of states to be stored in the variable
		 * @param i_dimensions The number of components in the variable
		 */
		variable (grids::grid &i_grid_m, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n ();
			data.resize (total * dimensions * i_states, 0.0);
			ld = 0;
			state = 0;
			last_update = 0;
			states = i_states;
		}

		/**
		 * @param i_grid_n A reference to the horizontal grid
		 * @param i_grid_m A reference to the vertical grid
		 * @param i_element_flags A reference to the global element flags
		 * @param i_name The string name of the variable
		 * @param i_states The number of states to be stored in the variable
		 * @param i_dimensions The number of components in the variable
		 */
		variable (grids::grid &i_grid_n, grids::grid &i_grid_m, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			grids.push_back (&i_grid_n);
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n () * i_grid_n.get_ld ();
			data.resize (total * i_states * dimensions, 0.0);
			ld = i_grid_m.get_n ();
			state = 0;
			last_update = 0;
			states = i_states;
		}

		/**
		 * @param n The number of dimensions
		 * @param i_grids An array of pointers to grid objects for the variable
		 * @param i_element_flags A reference to the global element flags
		 * @param i_name The string name of the variable
		 * @param i_states The number of states to be stored in the variable
		 * @param i_dimensions The number of components in the variable
		 */
		variable (int n, grids::grid **i_grids, int &i_element_flags, std::string i_name = "", int i_states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags), name (i_name) {
			total = n >= 1 ? 1 : 0;
			for (int i = 0; i < n; ++i)
			{
				grids.push_back (i_grids [i]);
				total *= i_grids [i]->get_ld ();
			}
			data.resize (total * i_states * dimensions, 0.0);
			ld = i_grids [n - 1]->get_ld ();
			state = 0;
			last_update = 0;
			states = i_states;
		}

		virtual ~variable () {}

		/**
		 * @brief Return a pointer to the requested state of the variable
		 * 
		 * @param i_state The integer state being requested
		 * @return A pointer to the data of the requested state
		 */
		double *ptr (int i_state = 0) {
			if (i_state >= states) {
				ERROR ("State " << i_state << " not initialized.");
				throw 501;
			}
			return &data [i_state * total * dimensions];
		}

		/**
		 * @return The name of the variable, constructing one if this one is compound and unnamed
		 */
		std::string get_name ();

		/**
		 * @return The size of one component (not one state) of the variable
		 */
		const int size () const {
			return total * dimensions;
		}

		/**
		 * @brief Return the number 
		 * @details [long description]
		 * @return [description]
		 */
		const int shape () const {
			return grids.size ();
		}

		/**
		 * @param n The index of the grid to exract
		 * @return A reference to the nth grid in the variable
		 */
		grid &get_grid (int n) {
			return *grids [n];
		}

		/**
		 * @return An array of pointers to the grids contained in the variable
		 */
		grid **get_grids () {
			return &grids [0];
		}

		/**
		 * @return The number of components in the variable
		 */
		const int dims () {
			return dimensions;
		}

		/**
		 * @return The number of states in the variable
		 */
		const int get_states () {
			return states;
		}

		/**
		 * @return The leading dimension of the data in the vertical dimension
		 */
		int get_ld () {
			return ld;
		}

		/**
		 * @brief Add a variable to this one such that it is constructed by arithmetic operations on other variables
		 * @details Since variables can be constructed by a series of arithmetic operations on other variables, additional variables can be added into this one. This is the method that should be called by operations on variable objects.
		 * 
		 * @param data A reference to the variable to add
		 * @param op The integer representation of the operation to engage for this variable
		 */
		void add_var (variable &data, int op) {
			inner.push_back (&data);
			vars.push_back (data.ptr ());
			ops.push_back (op);
			inner_states.push_back (-1);
		}

		/**
		 * @brief Reset the contents of any compound variables
		 * @details Set all the inner variable arrays to be empty.
		 */
		void reset_vars () {
			inner.resize (0);
			vars.resize (0);
			ops.resize (0);
			inner_states.resize (0);
		}

		/**
		 * @brief Update the value of the variable if it has any variable contents
		 * @details This will recursively call update on any contained variables as well.
		 * @return The new version of the variable
		 */
		int update ();

		/**
		 * @brief Index the 0th state of the variable
		 * @details This allows direct access to the variable contents and even marks the variable for updating; however, this is much slower than other methods of changing the variable contents.
		 * 
		 * @param index The integer index to index
		 */
		double &operator[] (int index) {
			component_flags &= ~updated;
			last_update = 0;
			state++;
			return data [index];
		}

		/**
		 * @brief Add another variable to the variable
		 * @details Generates a new compound variable that will update when either of its components change
		 * 
		 * @param other The variable to add to this variable
		 * @return A reference to a new compound variable
		 */
		variable &operator+ (variable &other) {
			std::shared_ptr <variable> new_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			new_var->add_var (*this, add);
			new_var->add_var (other, add);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Add a constant to the variable
		 * @details Generates two new variables, one a uniform variable and one a compound variable 
		 * 
		 * @param other The scalar to add to this variable
		 * @return A reference to a new compound variable
		 */
		variable &operator+ (double other) {
			std::shared_ptr <variable> new_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			std::shared_ptr <variable> uni_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, std::to_string ((long double) other), 1)));
			linalg::copy (this->size (), &other, uni_var->ptr (), 0);
			new_var->add_var (*this, add);
			new_var->add_var (*uni_var, add);
			store_var (new_var);
			store_var (uni_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Subtract another variable from the variable
		 * @details Generates a new compound variable that will update when either of its components change
		 * 
		 * @param other The variable to subtract from this variable
		 * @return A reference to a new compound variable
		 */
		variable &operator- (variable &other) {
			std::shared_ptr <variable> new_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			new_var->add_var (*this, add);
			new_var->add_var (other, sub);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Multiply another variable by the variable
		 * @details Generates a new compound variable that will update when either of its components change
		 * 
		 * @param other The variable to multiply by this variable
		 * @return A reference to a new compound variable
		 */
		variable &operator* (variable &other) {
			std::shared_ptr <variable> new_var;
			if (this->shape () > other.shape ()) {
				new_var.reset (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1));
			} else {
				new_var.reset (new variable (other.shape (), other.get_grids (), this->element_flags, "", 1));
			}			new_var->add_var (*this, add);
			new_var->add_var (other, mul);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Multiply the variable by a scalar
		 * @details Generates two new variables, one a uniform variable and one a compound variable 
		 * 
		 * @param other The scalar to multiply by this variable
		 * @return A reference to the new compound variable
		 */
		variable &operator* (double other) {
			std::shared_ptr <variable> new_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1)));
			std::shared_ptr <variable> uni_var (std::shared_ptr <variable> (new variable (this->shape (), this->get_grids (), this->element_flags, std::to_string ((long double) other), 1)));
			linalg::copy (this->size (), &other, uni_var->ptr (), 0);
			new_var->add_var (*this, add);
			new_var->add_var (*uni_var, mul);
			store_var (new_var);
			store_var (uni_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Divide this variable by another
		 * @details Generates a new compound variable that will update when either of its components change
		 * 
		 * @param other The variable to divide
		 * @return A reference to a new compound variable
		 */
		variable &operator/ (variable &other) {
			std::shared_ptr <variable> new_var;
			if (this->shape () > other.shape ()) {
				new_var.reset (new variable (this->shape (), this->get_grids (), this->element_flags, "", 1));
			} else {
				new_var.reset (new variable (other.shape (), other.get_grids (), this->element_flags, "", 1));
			}
			new_var->add_var (*this, add);
			new_var->add_var (other, div);
			store_var (new_var);
			new_var->update ();

			return *new_var;
		}

		/**
		 * @brief Reset the variable and set it to equal the right hand side
		 * 
		 * @param other The variable to set this one to equal
		 * @return A reference to this variable
		 */
		variable &operator== (variable &other) {
			this->reset_vars ();
			this->add_var (other, add);
			this->update ();
			return *this;
		}
	};

	/**
	 * @brief If the variable is the second argument in addition, make it the first
	 * 
	 * @param first A scalar
	 * @param other The variable to add this to
	 * @return A reference to a variable
	 */
	variable &operator+ (double first, variable &other);
	
	/**
	 * @brief If the variable is the second argument in subtraction, multiply by -1 then add it to the other
	 * 
	 * @param first A scalar
	 * @param other The variable to subtract this to
	 * @return A reference to a variable
	 */
	variable &operator- (double first, variable &other);

	/**
	 * @brief If the variable is the second argument in division and the first isn't, rework it to use the built-in operators
	 * 
	 * @param first A scalar
	 * @param other The variable to divide
	 * @return A reference to a variable
	 */
	variable &operator/ (double first, variable &other);
}

#endif /* end of include guard: VARIABLE_H__ */
