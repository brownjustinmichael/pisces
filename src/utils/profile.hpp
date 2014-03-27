/*!**********************************************************************
 * \file profile.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-03-26.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef PROFILE_HPP_J7WU4I34
#define PROFILE_HPP_J7WU4I34

#include "../config.hpp"

namespace utils
{
	enum profile_modes {
		average_profile = 0x01,
		median_profile = 0x02
	};
	
	template <class datatype>
	void profile (int m, int n, datatype *in, datatype *out, int mode = average_profile, int ldm = 0, int incout = 1) {
		if (!ldm) {
			ldm = m;
		}
		for (int i = 0; i < m; ++i) {
			if (mode == average_profile) {
				out [i * incout] = 0;
				for (int j = 0; j < m; ++j) {
					out [i * incout] += in [ldm * j + i];
				}
			} else {
				ERROR ("Unknown profiling mode");
				throw 0;
			}
		}
	}
} /* utils */


#endif /* end of include guard: PROFILE_HPP_J7WU4I34 */
