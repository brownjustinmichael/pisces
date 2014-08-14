/*!**********************************************************************
 * \file boundary.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-07-14.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef BOUNDARY_HPP_013D6464
#define BOUNDARY_HPP_013D6464

namespace bases
{
	template <class datatype>
	class boundary
	{
	public:
		boundary () {}
		
		virtual ~boundary () {}
		
		virtual int get_overlap () {
			return 0;
		}

		virtual int get_ex_overlap () {
			return 0;
		}
		
		virtual int get_excess () {
			return 0;
		}
		
		/*
			TODO Remove dimensionality from base class
		*/
		
		virtual void send (datatype *data_temp, int lda) {}
		
		virtual void receive (datatype *data_temp, int lda) {}
		
		virtual void calculate_rhs (datatype *data_in, datatype *interpolate_original, datatype *interpolate_data, datatype *data_out, int m, int lda, int flag) = 0;
		
		virtual void calculate_matrix (datatype timestep, datatype *default_matrix, datatype *matrix_in, datatype *interpolate_matrix, datatype *matrix_out, int lda, bool diverging = false) = 0;
	};
} /* bases */

#endif /* end of include guard: BOUNDARY_HPP_013D6464 */

