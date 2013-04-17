/*!***********************************************************************
 * \file boundary.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-09.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_CP5DR4CP
#define BOUNDARY_HPP_CP5DR4CP

#include "../config.hpp"
#include "../plan.hpp"

namespace boundary
{	
	class boundary : public plan
	{
	public:
		boundary (double i_alpha_plus, double *i_data_plus, double i_alpha_minus, double *i_data_minus) {
			alpha_plus = i_alpha_plus;
			alpha_minus = i_alpha_minus;
			data_plus = i_data_plus;
			data_minus = i_data_minus;
		}
		virtual ~boundary () {}
		
		virtual void execute () = 0;
	
	protected:
		double alpha_plus;
		double alpha_minus;
		double *data_plus;
		double *data_minus;
	};
	
	class boundary_1D : public boundary
	{
	public:
		boundary_1D (double i_alpha_plus, double *i_data_plus, double i_alpha_minus = 0, double *i_data_minus = NULL) : boundary (i_alpha_plus, i_data_plus, i_alpha_minus, i_data_minus) {}
		virtual ~boundary_1D () {}
		
		inline virtual void execute () {
			TRACE ("Executing...");
			
			DEBUG ("top = " << data_plus);
			DEBUG ("bottom = " << data_minus);
			
			if (!data_minus) {
				*data_plus *= alpha_plus;
			} else if (!data_plus) {
				*data_minus *= alpha_minus;
			} else {
				DEBUG ("before = " << *data_plus);
				*data_plus = alpha_plus * *data_plus + alpha_minus * *data_minus;
				DEBUG ("after = " << *data_plus);
				*data_minus = *data_plus;
			}
			
			TRACE ("Executed.")
		}
	};
} /* boundary */

#endif /* end of include guard: BOUNDARY_HPP_CP5DR4CP */
