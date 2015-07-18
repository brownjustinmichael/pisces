/*!**********************************************************************
 * \file real_plan.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-09-29.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef REAL_PLAN_HPP_CBB7844E
#define REAL_PLAN_HPP_CBB7844E

#include "grids/grid.hpp"
#include "linalg/utils.hpp"
#include "plan.hpp"

namespace plans
{
	/*!*******************************************************************
	 * \brief A subclass of plan, specific to real methods
	 * 
	 * These plans take input and produce output. Unlike implicit plans, they do not make any changes to the matrices associated with the solve. However, they do require a transform before the solve to put them in the correct state. This transform is not done within the class but rather in equation right before the solve happens (this is to make sure that the transform is done only once).
	 *********************************************************************/
	template <class datatype>
	class real_plan : public plan <datatype>
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
		real_plan (grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
		plans::plan <datatype> (i_data_in.ptr (), i_data_out.ptr (), &(i_data_in.element_flags), &(i_data_in.component_flags), i_coeff), 
		n (i_data_in.get_grid (0).get_n ()), 
		ldn (i_data_in.get_grid (0).get_ld ()), 
		m (i_data_in.get_grid (1).get_n ()), 
		dims (i_data_in.dims ()), 
		grid_n (i_data_in.get_grid (0)), 
		grid_m (i_data_in.get_grid (1)) {}

		virtual ~real_plan () {}
		
		virtual void setup () {}

		virtual int type () {
			return plan <datatype>::factory::real;
		}
		
		/*!*******************************************************************
		 * \copydoc plans::plan::execute ()
		 *********************************************************************/
		virtual void execute () = 0;
		
		/*!**********************************************************************
		 * \brief An abstract factory class designed to produce a real_plan instance
		 * 
		 * This factory class is implemented to make the creation of plan objects more intuitive. A factory can be given to an equation object, and the equation will know about the grids, data, and flags, so these need not be specified by the user.
		 ************************************************************************/
		class factory : public plan <datatype>::factory
		{
		public:
			factory (datatype i_coeff = 1.0) : plan <datatype>::factory (i_coeff) {}

			virtual ~factory () {}

			virtual const int type () const {
				return plan <datatype>::factory::real;
			}
		};
	};

	template <class datatype>
	class compound_plan : public real_plan <datatype>
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

		std::shared_ptr <grids::variable <datatype>> total, tmp;

	public:
		enum name
		{
			mult = 0x01,
			div = 0x02
		};

		compound_plan (grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out, datatype i_coeff = 1.0) : 
		real_plan <datatype> (i_data_in, i_data_out, i_coeff),
		data (i_data_in) {
			tmp = std::shared_ptr <grids::variable <datatype>> (new grids::variable <datatype> (i_data_in.get_grid (0), i_data_in.get_grid (1), *element_flags));
			total = std::shared_ptr <grids::variable <datatype>> (new grids::variable <datatype> (i_data_in.get_grid (0), i_data_in.get_grid (1), *element_flags));

		}

		~compound_plan () {}

		void add_plan (const typename plans::plan <datatype>::factory &i_factory, const int op = mult) {
			TRACE ("Adding plan...");
			if (i_factory.type () == plans::plan <datatype>::factory::impl) {
				throw 100;
			}
			plans.push_back (i_factory.instance (NULL, data, *tmp));
			operators.push_back (op);
		}

		void add_plan (std::shared_ptr <const typename plans::plan <datatype>::factory> i_factory, const int op = mult) {
			TRACE ("Adding plan...");
			add_plan (*i_factory, op);
		}
		
		void execute () {
			bool first = true;
			linalg::scale ((int) total_vec.size (), 0., total->ptr ());
			for (int i = 0; i < (int) plans.size (); ++i)
			{
				linalg::scale (tmp->size (), 0., tmp->ptr ());
				plans [i]->execute ();
				if (first) {
					linalg::add_scaled (total->size (), tmp->ptr (), total->ptr ());
					first = false;
					continue;
				}
				datatype *tmp_ptr = tmp->ptr ();
				datatype *total_ptr = tmp->ptr ();
				if (operators [i] == mult) {
					for (int j = 0; j < tmp->size (); ++j)
					{
						total_ptr [j] *= tmp_ptr [j];
					}
				} else if (operators [i] == div) {
					for (int j = 0; j < tmp->size (); ++j)
					{
						total_ptr [j] /= tmp_ptr [j];
					}
				} else {
					FATAL ("Unrecognized operator");
					throw 101;
				}
			}
			linalg::add_scaled (total->size (), total->ptr (), data_out);
		}

		class factory : public real_plan <datatype>::factory
		{
		protected:
			std::vector <std::shared_ptr <const typename plan <datatype>::factory>> factories;
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

			virtual std::shared_ptr <plan <datatype>> _instance (datatype **matrices, grids::variable <datatype> &i_data_in, grids::variable <datatype> &i_data_out) const {
				std::shared_ptr <compound_plan <datatype>> plan = std::shared_ptr <compound_plan <datatype>> (new compound_plan <datatype> (i_data_in, i_data_out));
				for (int i = 0; i < (int) factories.size (); ++i)
				{
					plan->add_plan (factories [i], operators [i]);
				}
				return plan;
			}
		};
	};

	std::shared_ptr <typename compound_plan <double>::factory> operator* (std::shared_ptr <typename plan <double>::factory> i_factory1, std::shared_ptr <typename plan <double>::factory> i_factory2);

	std::shared_ptr <typename compound_plan <double>::factory> operator/ (std::shared_ptr <typename plan <double>::factory> i_factory1, std::shared_ptr <typename plan <double>::factory> i_factory2);
} /* plans */

#endif /* end of include guard: REAL_PLAN_HPP_CBB7844E */
