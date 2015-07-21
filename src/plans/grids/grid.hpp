/*!***********************************************************************
 * \file grid.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-19.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef COLLOCATION_HPP_3FRTUK5Z
#define COLLOCATION_HPP_3FRTUK5Z

#include <vector>
#include <memory>
#include <cmath>

#include "versions/version.hpp"
#include "logger/logger.hpp"
#include "linalg/utils.hpp"

/*!**********************************************************************
 * \namespace grids
 * 
 * \brief A namespace containing grid objects with knowledge of the physical grid setup
 ************************************************************************/
namespace grids
{
	/*!**********************************************************************
	 * \brief A class containing the basic dimension information
	 ************************************************************************/
	class axis
	{
	private:
		int n; //!< The integer number of meshpoints along the axis
		int excess_0; //!< The integer number of excess meshpoints around index 0
		int excess_n; //!< The integer number of excess meshpoints around index n - 1
		double position_0; //!< The double position value at index excess_0
		double position_n; //!< The double position value at index n - 1 - excess_n
		
	public:
		/*!**********************************************************************
		 * \param i_n The integer number of meshpoints along the axis
		 * \param i_position_0 The double position value at index i_excess_0
		 * \param i_position_n The double position value at index n - 1 - i_excess_n
		 * \param i_excess_0 The integer number of excess meshpoints around index 0
		 * \param i_excess_n The integer number of excess meshpoints around index n - 1
		 ************************************************************************/
		axis (int i_n = 0, double i_position_0 = 0.0, double i_position_n = 0.0, int i_excess_0 = 0, int i_excess_n = 0) :
		n (i_n), excess_0 (i_excess_0), excess_n (i_excess_n), position_0 (i_position_0), position_n (i_position_n) {}
		
		/*!**********************************************************************
		 * \brief Get the total number of gridpoints along the axis
		 ************************************************************************/
		int get_n () {
			return n;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points near index 0
		 ************************************************************************/
		int get_excess_0 () {
			return excess_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points near index n - 1
		 ************************************************************************/
		int get_excess_n () {
			return excess_n;
		}
		
		/*!**********************************************************************
		 * \brief Get the position at index excess_0
		 ************************************************************************/
		double get_position_0 () {
			return position_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the position at index n - 1 - excess_n
		 ************************************************************************/
		double get_position_n () {
			return position_n;
		}
	};
	
	/*!*******************************************************************
	 * \brief A class containing a collocation grid
	 * 
	 * This will need to be instantiated for a collocation method. Ideally, the structure contains a grid of polynomial values at a number of collocation points which can be called for little temporal expense. Unfortunately, creating such a grid is an expensive procedure, so the creation is delayed until the get_data method is called.
	 *********************************************************************/
	template <class datatype>
	class grid
	{
	protected:
		int n; //!< The number of meshpoints along the axis
		int ld; //!< The number of meshpoints in the data along the axis (accounts for buffer points)
		int excess_0; //!< The number of excess meshpoints near index 0
		int excess_n; //!< The number of excess meshpoints near index n - 1
		datatype position_0; //!< The position at index excess_0
		datatype position_n; //!< The position at index n - 1 - excess_n
		int derivs; //!< The integer number of derivatives deep the collocation grid runs
		
		std::vector <datatype> positions; //!< A vector containing the positions along the axis
		bool calculated_matrix; //!< A boolean describing whether the collocation matrix has been calculated
	
	private:
		std::vector<std::vector<datatype> > data; //!< A double vector containing the vectors of collocation data

	public:
		/*!*******************************************************************
		 * \param i_axis_ptr A pointer to an axis object
		 * \param i_derivs The integer number of derivatives in the grid
		 * \param i_ld The extent of the data array
		 *********************************************************************/
		grid (axis *i_axis_ptr, int i_derivs, int i_ld = 0) :
		n (i_axis_ptr->get_n ()),
		ld (i_ld),
		excess_0 (i_axis_ptr->get_excess_0 ()),
		excess_n (i_axis_ptr->get_excess_n ()),
		position_0 (i_axis_ptr->get_position_0 ()),
		position_n (i_axis_ptr->get_position_n ()),
		derivs (i_derivs) {
			TRACE ("Instantiating...");
			if (ld == 0) {
				ld = n;
			}
			calculated_matrix = false;

			if (position_0 == position_n || n == 0) {
				ERROR ("Grid has no physical extent.");
				throw 0;
			}
			positions.resize (n + 1);

			TRACE ("Instantiated...");
		}
		
		/*!**********************************************************************
		 * \copydoc axis::axis
		 * 
		 * \param i_derivs The integer number of derivatives in the grid
		 * \param i_ld The extent of the data array
		 ************************************************************************/
		grid (int i_derivs, int i_n = 0, double i_position_0 = 0.0, double i_position_n = 0.0, int i_excess_0 = 0, int i_excess_n = 0, int i_ld = 0) {
			n = i_n;
			ld = i_ld;
			excess_0 = i_excess_0;
			excess_n = i_excess_n;
			position_0 = i_position_0;
			position_n = i_position_n;
			derivs = i_derivs;
			
			TRACE ("Instantiating...");
			if (ld == 0) {
				ld = n;
			}
			calculated_matrix = false;

			if (position_0 == position_n || n == 0) {
				ERROR ("Grid has no physical extent.");
				throw 0;
			}
			positions.resize (n + 1);

			TRACE ("Instantiated...");
		}
	
		virtual ~grid () {}
		
		/*!**********************************************************************
		 * \brief Return a reference to the position along the grid
		 * 
		 * \param index The integer index of interest
		 * 
		 * \return A reference to a position along the grid
		 ************************************************************************/
		const datatype& operator[] (int index) {
			return positions [index];
		}
		
		/*!**********************************************************************
		 * \brief Get the number of meshpoints along the axis
		 ************************************************************************/
		int get_n () {
			return n;
		}
		
		/*!**********************************************************************
		 * \brief Get the total array length along the axis (allowing for buffer points)
		 ************************************************************************/
		int get_ld () {
			return ld;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points by index 0
		 ************************************************************************/
		int get_excess_0 () {
			return excess_0;
		}
		
		/*!**********************************************************************
		 * \brief Get the number of overlapping points by index n - 1
		 ************************************************************************/
		int get_excess_n () {
			return excess_n;
		}
	
		/*!*******************************************************************
		 * \brief Get the data array for a given derivative level
		 * 
		 * \param deriv An integer derivative level (0 = value, 1 = first derivative, ...)
		 * 
		 * 
		 * This method is intended for use when a subroutine requires the matrix containing the collocation elements. If the matrix has not yet been calculated, this method will do so.
		 * 
		 * \return A pointer to the first element of the double data array
		 *********************************************************************/
		inline datatype* get_data (int deriv) {
			if (!calculated_matrix) {
				calculate_matrix ();
			}
			return &(data [deriv] [0]);
		}
		
	protected:
		/*!*******************************************************************
		 * \brief An indexing operation into the grid, for convenience
		 * 
		 * \param deriv The integer deriv to be indexed
		 * \param row The integer row to be indexed
		 * \param col The integer column to be indexed
		 * 
		 * \return The value at the index
		 *********************************************************************/
		inline datatype& index (int deriv, int row, int col) {
			return data [deriv] [row * ld + col];
		}
	
		/*!**********************************************************************
		 * \brief Calculate the collocation grid
		 * 
		 * This must be overwritten in subclasses that implement a grid.
		 ************************************************************************/
		virtual void _calculate_matrix () = 0;
		
	private:
		/*!**********************************************************************
		 * \brief Instantiate the collocation grid
		 * 
		 * Instantiate the collocation grid to be used by a collocation method. This method depends on the _calculate_matrix () operation.
		 ************************************************************************/
		virtual void calculate_matrix () {
			data.resize (derivs);
			for (int i = 0; i < derivs; ++i) {
				data [i].resize (ld * ld);
			}
			_calculate_matrix ();
			calculated_matrix = true;
		}
	};
	
#ifndef _VCOS
	/*!**********************************************************************
	 * \brief The mode of the simulation (vertical Chebyshev), to be used in output
	 ************************************************************************/
	enum mode {
		mode_flag = 0x10
	};
	
	/*!**********************************************************************
	 * \namespace grids::vertical
	 * 
	 * \brief A namespace containing grids in the vertical dimension
	 ************************************************************************/
	namespace vertical
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 * 
		 * This collocation grid stores the N collocation points for up to the Mth order Chebyshev polynomial and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public grids::grid <datatype>
		{
		private:
			using grids::grid <datatype>::n;
			using grids::grid <datatype>::ld;
			using grids::grid <datatype>::excess_0;
			using grids::grid <datatype>::excess_n;
			using grids::grid <datatype>::position_0;
			using grids::grid <datatype>::position_n;
			using grids::grid <datatype>::derivs;
			using grids::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			std::vector<bool> exists_array; //!< A bool vector containing whether the points exist
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
			
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
			
			/*!**********************************************************************
			 * \copydoc axis::axis
			 * 
			 * \param i_ld The extent of the data array
			 ************************************************************************/
			grid (int i_n = 0, double i_position_0 = 0.0, double i_position_n = 0.0, int i_excess_0 = 0, int i_excess_n = 0, int i_ld = 0);
			
			virtual ~grid () {};
			
		protected:
			/*!**********************************************************************
			 * \copydoc grids::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
			
		private:
			/*!*******************************************************************
			* \brief A check to see whether an index exists, for recursion
			* 
			* \param d The integer deriv to be indexed
			* \param m The integer row to be indexed
			* \param k The integer column to be indexed
			* 
			* \return An std::vector<bool>::reference of whether the element exists
			*********************************************************************/
			std::vector<bool>::reference exists (int d, int m, int k) {
				return exists_array [d + m * derivs + k * derivs * n];
			}
	
			/*!*******************************************************************
			 * \brief A recursive method to calculate the full collocation grid
			 * 
			 * \param d The integer deriv to be indexed
			 * \param m The integer row to be indexed
			 * \param k The integer column to be indexed
			 * 
			 * Uses known relationships of Chebyshev polynomials to calculate the value of the polynomial and its derivatives at an arbitrary collocation grid point.
			 * 
			 * \return The datatype value of the index
			 *********************************************************************/
			datatype recursion (int d, int m, int k);
		};
	} /* vertical */
#else
	/*!**********************************************************************
	 * \brief The mode of the simulation (vertical cosine), to be used in output
	 ************************************************************************/
	enum mode {
		mode_flag = 0x20
	};
	
	
	// A note: The cosine expansion doesn't work because it restricts the derivative at the inner boundaries, this would need to be fourier
	namespace vertical
	{
		/*!*******************************************************************
		 * \brief A collocation grid for cosines
		 *
		 * This collocation grid stores the N collocation points for up to the Mth order cosine modes and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public grids::grid <datatype>
		{
		private:
			using grids::grid <datatype>::n;
			using grids::grid <datatype>::ld;
			using grids::grid <datatype>::excess_0;
			using grids::grid <datatype>::excess_n;
			using grids::grid <datatype>::position_0;
			using grids::grid <datatype>::position_n;
			using grids::grid <datatype>::derivs;
			using grids::grid <datatype>::positions;

			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations

		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
			
			virtual ~grid () {};

		protected:
			/*!**********************************************************************
			 * \copydoc grids::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
		};
	} /* vertical */
#endif
	
	/*!**********************************************************************
	 * \namespace grids::horizontal
	 * 
	 * \brief A namespace containing grids in the horizontal dimension
	 ************************************************************************/
	namespace horizontal
	{
		/*!*******************************************************************
		 * \brief A collocation grid for Chebyshev polynomials
		 * 
		 * This collocation grid stores the N collocation points for up to the Mth order fourier mode and their first and second derivatives
		 *********************************************************************/
		template <class datatype>
		class grid : public grids::grid <datatype>
		{
		private:
			using grids::grid <datatype>::n;
			using grids::grid <datatype>::ld;
			using grids::grid <datatype>::excess_0;
			using grids::grid <datatype>::excess_n;
			using grids::grid <datatype>::position_0;
			using grids::grid <datatype>::position_n;
			using grids::grid <datatype>::derivs;
			using grids::grid <datatype>::positions;
	
			datatype scale; //!< A datatype by which the collocation grid should be scaled
			datatype width; //!< The datatype width of the collocation region
			datatype pioN; //!< The datatype 3.14159.../N, for use in calculations
			
		public:
			/*!*******************************************************************
			 * \param i_axis_ptr A pointer to an axis object
			 *********************************************************************/
			grid (axis *i_axis_ptr);
			
			/*!**********************************************************************
			 * \copydoc axis::axis
			 * 
			 * \param i_ld The extent of the data array
			 ************************************************************************/
			grid (int i_n = 0, double i_position_0 = 0.0, double i_position_n = 0.0, int i_excess_0 = 0, int i_excess_n = 0, int i_ld = 0);
				
			virtual ~grid () {};
			
			/*!**********************************************************************
			 * \return The version of the class
			 ************************************************************************/
			static versions::version& version () {
				static versions::version version ("1.1.0.0");
				return version;
			}
			
		protected:
			/*!**********************************************************************
			 * \copydoc grids::grid::_calculate_matrix()
			 ************************************************************************/
			void _calculate_matrix ();
		};
	} /* fourier */

	template <class datatype>
	class compound_variable;

	template <class datatype>
	class variable
	{
	protected:
		static std::vector <std::shared_ptr <variable <datatype>>> tmps;

		std::vector <grids::grid <datatype> *> grids;
		int dimensions;
		std::vector <datatype> data;
		std::vector <datatype *> vars;
		// TODO This won't work. We need a way to keep references or entire variables
		std::vector <variable <datatype> *> inner;
		std::vector <int> ops;
		bool needs_update = false;
		int ld;
		int total;

	public:
		int last_update;
		int component_flags;
		int &element_flags;

		enum operators {
			add = 0x01,
			sub = 0x02,
			mul = 0x03,
			div = 0x04
		};

		variable (grids::grid <datatype> &i_grid_m, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			grids.push_back (&i_grid_m);
			total = i_grid_m.get_n ();
			data.resize (total * dimensions * states, 0.0);
			ld = 0;
		}

		variable (grids::grid <datatype> &i_grid_n, grids::grid <datatype> &i_grid_m, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			grids.push_back (&i_grid_n);
			grids.push_back (&i_grid_m);
			data.resize (i_grid_m.get_n () * i_grid_n.get_ld () * states * dimensions, 0.0);
			ld = i_grid_m.get_n ();
			total = i_grid_m.get_n () * i_grid_n.get_ld ();
		}

		variable (int n, grids::grid <datatype> **i_grids, int &i_element_flags, int states = 3, int i_dimensions = 1) : dimensions (i_dimensions), component_flags (0x00), element_flags (i_element_flags) {
			total = n > 1 ? 1 : 0;
			for (int i = 0; i < n; ++i)
			{
				grids.push_back (i_grids [i]);
				total *= i_grids [i]->get_ld ();
			}
			data.resize (total * states * dimensions, 0.0);
			ld = i_grids [n - 1]->get_ld ();
		}

		virtual ~variable () {}

		datatype *ptr (int state = 0) {
			return &data [state * total * dimensions];
		}

		const int size () const {
			return total * dimensions;
		}

		const int shape () const {
			return grids.size ();
		}

		grid <datatype> &get_grid (int n) {
			return *grids [n];
		}

		grid <datatype>** get_grids () {
			return &grids [0];
		}

		const int dims () {
			return dimensions;
		}

		const int get_ld () {
			return ld;
		}

		void add_var (variable <datatype> &data, int op) {
			needs_update = true;
			DEBUG ("ADDING " << data.ptr ());
			inner.push_back (&data);
			vars.push_back (data.ptr ());
			ops.push_back (op);
		}

		void reset_vars () {
			needs_update = false;
			inner.resize (0);
			vars.resize (0);
			ops.resize (0);
		}

		bool update ();

		variable <datatype> &operator== (variable <datatype> &other) {
			this->reset_vars ();
			this->add_var (other, add);
			return *this;
		}

		variable <datatype> &operator* (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags)));
			new_var->add_var (*this, add);
			new_var->add_var (other, mul);
			tmps.push_back (new_var);

			return *new_var;
		}

		variable <datatype> &operator/ (variable <datatype> &other) {
			std::shared_ptr <variable <datatype>> new_var (std::shared_ptr <variable <datatype>> (new variable <datatype> (this->shape (), this->get_grids (), this->element_flags)));
			new_var->add_var (*this, add);
			new_var->add_var (other, div);
			tmps.push_back (new_var);

			return *new_var;
		}
	};
} /* grids */

#endif /* end of include guard: COLLOCATION_HPP_3FRTUK5Z */
