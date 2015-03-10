/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>
#include "io/input.hpp"
#include "io/output.hpp"
#include "io/formats/netcdf.hpp"
#include "io/formats/virtual.hpp"

#define TEST_TINY 1.e-4

class io_two_d_test_suite : public CxxTest::TestSuite
{
public:
	void test_formatted_netcdf () {
		int n = 100, m = 200;
		std::vector <double> init (n * m), init_copy (n * m), profile (n), profile2 (m);
		double scalar, scalar_copy;
		std::string file_name = "output";
		{
			io::formatted_output <formats::netcdf> output_stream (formats::data_grid::two_d (n, m), file_name, formats::replace_file);
			
			srand (1);
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					init [i * m + j] = rand () % 100;
					init_copy [i * m + j] = init [i * m + j];
				}
			}
			scalar = rand () % 100;
			scalar_copy = scalar;
	
			output_stream.append <double> ("test", &init [0]);
			output_stream.append <double> ("profile", &init [0], formats::one_d);
			output_stream.append <double> ("profile2", &init [0], formats::m_profile);
			output_stream.append <double> ("scale", &scalar, formats::scalar);
			output_stream.to_file ();
		}

		for (int i = 0; i < n * m; ++i) {
			init [i] = 0.0;
		}
		scalar = 0.0;
	
		io::formatted_input <formats::netcdf> input_stream (formats::data_grid::two_d (n, m), file_name);
	
		input_stream.append <double> ("test", &init [0]);
		input_stream.append <double> ("profile2", &profile2 [0], formats::m_profile);
		input_stream.append <double> ("profile", &profile [0], formats::one_d);
		input_stream.append <double> ("scale", &scalar, formats::scalar);
		input_stream.from_file ();
	
		for (int i = 0; i < n; ++i) {
			TSM_ASSERT_DELTA ("IO failure in formatted netcdf", profile [i], init_copy [i], TEST_TINY);
		}

		for (int i = 0; i < m; ++i) {
			TSM_ASSERT_DELTA ("IO failure in formatted netcdf", profile2 [i], init_copy [i], TEST_TINY);
		}

		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("IO failure in formatted netcdf", init [i * m + j], init_copy [i * m + j], TEST_TINY);
			}
		}
		TSM_ASSERT_DELTA ("IO scalar failure in formatted netcdf", scalar, scalar_copy, TEST_TINY);
		std::cout << "COMPLETE";
	}
	
	void test_appending_netcdf () {
		int n = 100, m = 200, records = 10;
	
		std::string file_name = "output_appending";
	
		std::vector <double> init (n * m * records), init_copy (n * m * records);
		std::vector <double> scalar (records), scalar_copy (records);
	
		srand (1);
		{
			formats::data_grid io_grid = formats::data_grid::two_d (n, m);
			std::cout << io_grid.get_n_dims ();
			io::appender_output <formats::netcdf> output_stream (formats::data_grid::two_d (n, m), file_name, 1);
			for (int k = 0; k < records; ++k) {
				output_stream.append <double> ("test", &init [k * n * m]);
				output_stream.append <double> ("scale", &scalar [k], formats::scalar);
				for (int i = 0; i < n; ++i) {
					for (int j = 0; j < m; ++j) {
						init [k * n * m + i * m + j] = rand () % 100;
						init_copy [k * n * m + i * m + j] = init [k * n * m + i * m + j];
					}
				}
				scalar [k] = rand () % 100;
				scalar_copy [k] = scalar [k];
				
				output_stream.to_file ();
			
				for (int i = 0; i < n * m; ++i) {
					init [i] = 0.0;
				}
				scalar [k] = 0.0;
			}
		}
		
		io::formatted_input <formats::netcdf> input_stream (formats::data_grid::two_d (n, m), file_name);
		
		for (int k = 0; k < records; ++k) {
			input_stream.append <double> ("test", &init [k * n * m]);
			input_stream.append <double> ("scale", &scalar [k], formats::scalar);
			input_stream.from_file (k);
				
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					TSM_ASSERT_DELTA ("IO failure in appender netcdf", init [k * n * m + i * m + j], init_copy [k * n * m + i * m + j], TEST_TINY);
				}
			}
			TSM_ASSERT_DELTA ("IO scalar failure in appender netcdf", scalar [k], scalar_copy [k], TEST_TINY);
			
		}
	}
	
	void test_virtual_file () {
		int n = 100, m = 200;
	
		std::string file_name = "output_virtual";
	
		std::vector <double> init (n * m), init_copy (n * m);
	
		srand (1);
	
		io::appender_output <formats::virtual_format> output_stream (formats::data_grid::two_d (n, m), file_name, 1);
		output_stream.append <double> ("test", &init [0]);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				init [i * m + j] = rand () % 100;
				init_copy [i * m + j] = init [i * m + j];
			}
		}
		
		output_stream.to_file ();
		
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("IO failure in virtual format output", formats::virtual_files [file_name].index <double> ("test", i, j), init_copy [i * m + j], TEST_TINY);
			}
		}
	
		for (int i = 0; i < n * m; ++i) {
			init [i] = 0.0;
		}
		
		io::formatted_input <formats::virtual_format> input_stream (formats::data_grid::two_d (n, m), file_name);
		
		input_stream.append <double> ("test", &init [0]);
		input_stream.from_file ();
			
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("IO failure in virtual format input", init [i * m + j], init_copy [i * m + j], TEST_TINY);
			}
		}
	}
};
