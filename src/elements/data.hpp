/*!**********************************************************************
 * \file data.hpp
 * /Users/justinbrown/Dropbox/pisces
 * 
 * Created by Justin Brown on 2014-10-09.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef DATA_HPP_574F2B9E
#define DATA_HPP_574F2B9E

#include <string>
#include <map>
#include <memory>

#include "linalg/utils.hpp"

#include "io/input.hpp"
#include "io/output.hpp"
#include "io/formats/exceptions.hpp"
#include "plans/grid.hpp"
#include "plans-transforms/implemented_transformer.hpp"

namespace data
{
	enum initialize_element_flags {
		uniform_n = 0x01,
		uniform_m = 0x02
	};
	
	template <class datatype>
	class data
	{
	private:
		std::map <int, std::string> scalar_names;
		static int mode;
		
		std::shared_ptr <io::output> normal_stream; //!< An implementation to output in normal space
		std::shared_ptr <io::output> transform_stream; //!< An implementation to output in transform space
		std::shared_ptr <io::output> stat_stream; //!< An implementation to output in transform space
		std::vector <std::shared_ptr <io::output>> streams;
		std::vector <int> stream_conditions;
		std::vector <bool> done;
		std::vector <std::shared_ptr <io::output>> normal_profiles; //!< An implementation to output in transform space
		
	public:
		datatype duration;
		datatype timestep;
		std::map <int, int> flags;
		data (int i_transform_threads = 1) : transform_threads (i_transform_threads) {
			flags [0] = 0x0;
			duration = 0.0;
		}
		
		virtual ~data () {}
		
		datatype *operator[] (const int name) {
			return &(scalars [name] [0]);
		}
		
		datatype *operator() (const int name, const int index = 0) {
			return &(scalars [name] [index]);
		}
		
		const int &get_mode () {
			return mode;
		}
		
		void output () {
			INFO ("OUTPUT");
			if (normal_stream) normal_stream->to_file ();
		}
		
		void output_transform () {
			if (transform_stream) transform_stream->to_file ();
		}
		
		void output_stat () {
			if (stat_stream) stat_stream->to_file ();
		}
		
		virtual datatype *initialize (int i_name, std::string i_str, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			flags [i_name] = 0x00;
			scalar_names [i_name] = i_str;
			return _initialize (i_name, initial_conditions, i_flags);
		}
		
		void transform (int i_flags) {
			TRACE ("Transforming...");
			int threads = transform_threads;
			for (int i = 0; i < (int) streams.size (); ++i) {
				if (done [i]) continue;
				if (stream_conditions [i] & transformed_vertical) {
					if (!(flags [transforms [0]] & transformed_vertical)) {
						continue;
					}
				} else {
					if ((flags [transforms [0]] & transformed_vertical)) {
						continue;
					}
				}
				if (stream_conditions [i] & transformed_horizontal) {
					if (!(flags [transforms [0]] & transformed_horizontal)) {
						continue;
					}
				} else {
					if ((flags [transforms [0]] & transformed_horizontal)) {
						continue;
					}
				}
				streams [i]->to_file ();
				done [i] = true;
			}
#pragma omp parallel num_threads (threads)
			{
#pragma omp for
				for (int i = 0; i < (int) transforms.size (); ++i) {
					transformers [transforms [i]]->transform (i_flags);
				}
			}
		}
		
		void reset () {
			for (int i = 0; i < (int) done.size (); ++i) {
				done [i] = false;
			}
		}
		
		/*!**********************************************************************
		 * \brief Given an input stream, load all the relevant data into the element
		 * 
		 * \param input_stream_ptr A pointer to an input object
		 ************************************************************************/
		void setup (io::input *input_ptr) {
			// Iterate through the scalar fields and append them to the variables for which the input will search
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				input_ptr->template append <datatype> (scalar_names [iter->first], (*this) (iter->first));
			}
			
			input_ptr->template append <datatype> ("t", &duration, io::scalar);
			int mode;
			input_ptr->template append <int> ("mode", &mode, io::scalar);
			
			try {
				// Read from the input into the element variables
				input_ptr->from_file ();
			} catch (io::formats::exceptions::bad_variables &except) {
				// Send a warning if there are missing variables in the input
				WARN (except.what ());
			}

			// Make sure that the mode matched the mode of the element
			if (mode != get_mode ()) {
				FATAL ("Loading simulation in different mode: " << mode << " instead of " << get_mode ());
				throw 0;
			}
		}
	
		/*!**********************************************************************
		 * \brief Given an output_stream, append all relevant outputs
		 * 
		 * \param output_ptr A shared_ptr to the output object
		 * \param flags The integer flags to specify the type of output (transform stream, normal_stream)
		 * 
		 * Given an output object, prepare it to output all the scalar fields tracked directly by the element. Make it one of the output streams called during output. If flags is normal_stream, output the variables when in Cartesian space, else if flags is transform_stream, output with horizontal grid in Fourier space, but vertical grid in Cartesian space.
		 ************************************************************************/
		virtual void setup_output (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Iterate through the scalar fields and append them to the variables which the output will write to file
			typedef typename std::map <int, std::vector <datatype> >::iterator iterator;
			for (iterator iter = scalars.begin (); iter != scalars.end (); ++iter) {
				output_ptr->template append <datatype> (scalar_names [iter->first], (*this) (iter->first));
			}
			
			output_ptr->template append <datatype> ("t", &duration, io::scalar);
			output_ptr->template append <const int> ("mode", &(get_mode ()), io::scalar);

			// Check the desired output time and save the output object in the appropriate variable
			streams.push_back (output_ptr);
			done.push_back (false);
			// stream_conditions.push_back (flags);
			if (flags & transform_output) {
				stream_conditions.push_back (transformed_vertical | transformed_horizontal);
				// transform_stream = output_ptr;
			} else if (flags & normal_output) {
				stream_conditions.push_back (0x00);
				// normal_stream = output_ptr;
			}
		}
		
		virtual void setup_stat (std::shared_ptr <io::output> output_ptr, int flags = 0x00) {
			// Also prepare to output the total simulated time and geometry mode
			output_ptr->template append <datatype> ("t", &duration, io::scalar);
			output_ptr->template append <datatype> ("dt", &timestep, io::scalar);
			// Check the desired output time and save the output object in the appropriate variable
			DEBUG ("SETTING STAT")
			stat_stream = output_ptr;
		}
		
	protected:
		std::map <int, std::vector <datatype>> scalars;
		std::map <int, std::shared_ptr <plans::transformer <datatype>>> transformers;
		virtual datatype *_initialize (int i_name, datatype* initial_conditions = NULL, int i_flags = 0x00) = 0;
		std::vector <int> transforms; //!< A vector of integer keys to the transform maps
		int transform_threads;
	};
	
	template <class datatype>
	class implemented_data : public data <datatype>
	{
	private:
		std::shared_ptr <plans::grid <datatype>> grid_n, grid_m;
		int n, m;
		
	public:
		implemented_data (plans::axis *i_axis_n, plans::axis *i_axis_m, int i_transform_threads = 1) : data <datatype> (i_transform_threads), grid_n (std::shared_ptr <plans::grid <datatype>> (new typename plans::horizontal::grid <datatype> (i_axis_n))), grid_m (std::shared_ptr <plans::grid <datatype>> (new typename plans::vertical::grid <datatype> (i_axis_m))), n (grid_n->get_n ()), m (grid_m->get_n ()) {}
		
		virtual ~implemented_data () {}
		
		datatype *operator () (const int name, const int i = 0, const int j = 0) {
			return data <datatype>::operator() (name, i * grid_m->get_ld () + j); 
		}
	
	protected:
		virtual datatype *_initialize (int name, datatype* initial_conditions = NULL, int i_flags = 0x00) {
			TRACE ("Initializing " << name << "...");
			// Size allowing for real FFT buffer
			data <datatype>::scalars [name].resize (grid_n->get_ld () * m, 0.0);
			if (name == x_position) {
				for (int j = 0; j < m; ++j) {
					linalg::copy (n, &((*grid_n) [0]), (*this) (name, 0, j), 1, m);
				}
			} else if (name == z_position) {
				for (int i = 0; i < n; ++i) {
					linalg::copy (m, &((*grid_m) [0]), (*this) (name, i));
				}
			} else {
				if (initial_conditions) {
					if ((i_flags & uniform_m) && (i_flags & uniform_n)) {
						for (int i = 0; i < n; ++i) {
							for (int j = 0; j < m; ++j) {
								*(*this) (name, i, j) = *initial_conditions;
							}
						}
					} else if (i_flags & uniform_m) {
						for (int j = 0; j < m; ++j) {
							linalg::copy (n, initial_conditions, (*this) (name, 0, j), 1, m);
						}
					} else if (i_flags & uniform_n) {
						for (int i = 0; i < n; ++i) {
							linalg::copy (m, initial_conditions, (*this) (name, i));
						}
					} else {
						linalg::copy (n * m, initial_conditions, (*this) (name));
					}
				}
			}
			
			if ((name != x_position) && (name != z_position)) {
				data <datatype>::transforms.push_back (name);
				data <datatype>::transformers [name] = std::shared_ptr <plans::transformer <datatype> > (new plans::implemented_transformer <datatype> (*grid_n, *grid_m, (*this) (name), NULL, forward_vertical | forward_horizontal | inverse_vertical | inverse_horizontal , &data <datatype>::flags [state], &data <datatype>::flags [name], data <datatype>::transform_threads));
			}
			TRACE ("Done.");
			return &(data <datatype>::scalars [name] [0]);
		}
	};
} /* data */

#endif /* end of include guard: DATA_HPP_574F2B9E */
