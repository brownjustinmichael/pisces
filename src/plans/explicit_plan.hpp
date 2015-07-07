/*!**********************************************************************
 * \file explicit_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef EXPLICIT_PLAN_HPP_5171CAFD
#define EXPLICIT_PLAN_HPP_5171CAFD

#include "grids/grid.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to explicit methods
	 * 
	 * These plans take input and produce output. Unlike implicit plans, they do not alter the matrix that is solved. Their output also does not need to be retransformed before the solve happens, unlike real plans.
	 *********************************************************************/
	template <class datatype>
	class explicit_plan : public plan <datatype>
	{
	protected:
		using plan <datatype>::coeff;
		using plan <datatype>::data_in;
		using plan <datatype>::data_out;
		int n; //!< An integer number of data elements (grid points) in the horizontal
		int ldn; //!< The integer max dimension in the horizontal direction
		int m; //!< The integer number of data elements in the vertical
		int dims;
		grids::grid <datatype> &grid_n; //!< A reference to the horizontal grid object
		grids::grid <datatype> &grid_m; //!< A reference to the vertical grid object
		
	public:
		/*!**********************************************************************
		 * \param i_grid_n The grid object associated with the horizontal direction
		 * \param i_grid_m The grid object associated with the vertical direction
		 * \param i_data_in The input data for the plan
		 * \param i_data_out The output data for the plan
		 * \param i_element_flags A pointer to integer flags associated with the element on the whole
		 * \param i_component_flags A pointer to the integer flags associated with the variable associated with the plan
		 ************************************************************************/
		explicit_plan (grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, datatype i_coeff = 1.0, int *i_element_flags = NULL, int *i_component_flags = NULL) : 
		plans::plan <datatype> (i_data_in.ptr (), i_data_out, i_element_flags, i_component_flags, i_coeff), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}

		virtual ~explicit_plan () {}
		
		virtual void setup () {}
		
		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
	
		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce an explicit_plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory : public plan <datatype>::factory
		{
		public:
			factory (datatype i_coeff = 1.0) : plan <datatype>::factory (i_coeff) {}

			virtual ~factory () {}
			
			virtual const int type () const {
				return plan <datatype>::factory::expl;
			}
		};
	};

	template <class datatype>
	class compound_plan : public explicit_plan <datatype>
	{
	protected:
		using plan <datatype>::data_in;
		using plan <datatype>::data_out;
		using plan <datatype>::element_flags;
		using plan <datatype>::component_flags;

		grids::variable <datatype> &data;

		std::vector <std::shared_ptr <plan <datatype>>> plans;
		std::vector <int> operators;
		std::vector <datatype> total_vec;
		std::vector <datatype> tmp_vec;

		datatype *total, tmp;

	public:
		enum name
		{
			mult = 0x01,
			div = 0x02
		};

		compound_plan (grids::variable <datatype> *i_data_in, datatype *i_data_out, int *i_element_flags = NULL, int *i_component_flags = NULL, datatype i_coeff = 1.0) : 
		explicit_plan <datatype> (i_data_in, i_data_out, i_element_flags, i_component_flags, i_coeff),
		data (i_data_in) {
			tmp_vec.resize (i_data_in.size ());
			total_vec.resize (i_data_in.size ());
			tmp = &tmp_vec [0];
			total = &total_vec [0];
		}

		~compound_plan () {}

		void add_plan (const typename plans::plan <datatype>::factory &i_factory, int op = mult) {
			TRACE ("Adding plan...");
			if (i_factory.type () != plans::plan <datatype>::factory::expl) {
				throw 100;
			}
			plans.push_back (i_factory.instance (NULL, data, &tmp [0], element_flags, component_flags));
			operators.push_back (op);
		}
		
		void execute () {
			bool first = true;
			for (int i = 0; i < (int) plans.size (); ++i)
			{
				plans [i].execute ();
				if (first) {
					linalg::add_scaled ((int) total_vec.size (), tmp, total);
					first = false;
					continue;
				}
				if (operators [i] == mult) {
					for (int j = 0; j < (int) tmp_vec.size (); ++j)
					{
						total [j] *= tmp [j];
					}
				} else if (operators [i] == div) {
					for (int j = 0; j < (int) tmp_vec.size (); ++j)
					{
						total [j] /= tmp [j];
					}
				} else {
					FATAL ("Unrecognized operator");
					throw 101;
				}
			}
			linalg::add_scaled ((int) total_vec.size (), total, data_out);
		}

		class factory : public plan <datatype>::factory
		{
		protected:
			std::vector <std::shared_ptr <typename plan <datatype>::factory>> factories;
			std::vector <int> operators;

		public:
			factory () {}

			~factory () {}

			void add_plan (std::shared_ptr <typename plan <datatype>::factory> i_factory, int op = mult) {
				factories.push_back (i_factory);
				operators.push_back (op);
			}

			std::shared_ptr <typename plan <datatype>::factory> operator* (std::shared_ptr <typename plan <datatype>::factory> i_factory) {
				add_plan (i_factory, mult);
			}

			std::shared_ptr <typename plan <datatype>::factory> operator/ (std::shared_ptr <typename plan <datatype>::factory> i_factory) {
				add_plan (i_factory, div);
			}

			virtual std::shared_ptr <plan <datatype>> _instance (datatype **matrices, grids::variable <datatype> &i_data_in, datatype *i_data_out = NULL, int *i_element_flags = NULL, int *i_component_flags = NULL) const {
				std::shared_ptr <plan <datatype>> plan = std::shared_ptr <compound_plan <datatype>> (new compound_plan <datatype> (i_data_in, i_data_out, i_element_flags, i_component_flags));
				for (int i = 0; i < (int) factories.size (); ++i)
				{
					add_plan (factories [i], operators [i]);
				}
				return plan;
			}
		};
	};

	template <class datatype>
	std::shared_ptr <typename compound_plan <datatype>::factory> operator* (std::shared_ptr <typename plan <datatype>::factory> i_factory1, std::shared_ptr <typename plan <datatype>::factory> i_factory2) {
		std::shared_ptr <typename compound_plan <datatype>::factory> plan = typename std::shared_ptr <typename compound_plan <datatype>::factory> (new typename compound_plan <datatype>::factory ());

		plan->add_plan (i_factory1, compound_plan <datatype>::mult);
		plan->add_plan (i_factory2, compound_plan <datatype>::mult);
		return plan;
	}

	template <class datatype>
	std::shared_ptr <typename compound_plan <datatype>::factory> operator/ (std::shared_ptr <typename plan <datatype>::factory> i_factory1, std::shared_ptr <typename plan <datatype>::factory> i_factory2) {
		std::shared_ptr <typename compound_plan <datatype>::factory> plan = std::shared_ptr <typename compound_plan <datatype>::factory> (new typename compound_plan <datatype>::factory ());

		plan->add_plan (i_factory1, compound_plan <datatype>::mult);
		plan->add_plan (i_factory2, compound_plan <datatype>::div);
		return plan;
	}
} /* plans */

#endif /* end of include guard: EXPLICIT_PLAN_HPP_5171CAFD */
