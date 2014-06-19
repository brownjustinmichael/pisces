/*!**********************************************************************
 * \file transform_test.cpp
 * /Developer/NVIDIA/CUDA-5.5/samples/7_CUDALibraries/simpleCUFFT
 * 
 * Created by Justin Brown on 2014-03-01.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#include <cxxtest/TestSuite.h>
#include "../../src/utils/io.hpp"
#include "../../src/utils/formats.hpp"

#define TEST_TINY 1.e-4

class io_two_d_test_suite : public CxxTest::TestSuite
{
public:
	void test_formatted_netcdf () {
		int n = 100, m = 200;
		std::vector <double> init (n * m), init_copy (n * m);
		double scalar, scalar_copy;
		std::string file_name = "output";
		{
			io::formatted_output <io::formats::two_d::netcdf> output_stream (file_name, io::replace_file, n, m, 1, 0, 0, 0, 0, 0);
	
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
			output_stream.append_scalar <double> ("scale", &scalar);
			output_stream.to_file ();
		}

		for (int i = 0; i < n * m; ++i) {
			init [i] = 0.0;
		}
		scalar = 0.0;
	
		io::formatted_input <io::formats::two_d::netcdf> input_stream (file_name, n, m, 1, 0, 0, 0, 0, 0);
	
		input_stream.append <double> ("test", &init [0]);
		input_stream.append_scalar <double> ("scale", &scalar);
		input_stream.from_file ();
	
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("IO failure in formatted netcdf", init [i * m + j], init_copy [i * m + j], TEST_TINY);
			}
		}
		TSM_ASSERT_DELTA ("IO scalar failure in formatted netcdf", scalar, scalar_copy, TEST_TINY);
	}
	
	void test_appending_netcdf () {
		int n = 100, m = 200, records = 10;
	
		std::string file_name = "output_appending";
	
		std::vector <double> init (n * m * records), init_copy (n * m * records);
		std::vector <double> scalar (records), scalar_copy (records);
	
		srand (1);
	
		{
			io::appender_output <io::formats::two_d::netcdf> output_stream (file_name, 1, n, m, 1, 0, 0, 0, 0, 0);
			for (int k = 0; k < records; ++k) {
				output_stream.append <double> ("test", &init [k * n * m]);
				output_stream.append_scalar <double> ("scale", &scalar [k]);
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
		
		io::formatted_input <io::formats::two_d::netcdf> input_stream (file_name, n, m, 1, 0, 0, 0, 0, 0);
		
		for (int k = 0; k < records; ++k) {
			input_stream.append <double> ("test", &init [k * n * m]);
			input_stream.append_scalar <double> ("scale", &scalar [k]);
			input_stream.from_file (k);
				
			for (int i = 0; i < n; ++i) {
				for (int j = 0; j < m; ++j) {
					TSM_ASSERT_DELTA ("IO failure in appender netcdf", init [k * n * m + i * m + j], init_copy [k * n * m + i * m + j], TEST_TINY);
				}
			}
			TSM_ASSERT_DELTA ("IO scalar failure in appender netcdf", scalar [k], scalar_copy [k], TEST_TINY);
		}
	}
	
	void test_virtual_dump () {
		int n = 100, m = 200;
	
		std::string file_name = "output_virtual";
	
		std::vector <double> init (n * m), init_copy (n * m);
	
		srand (1);
	
		io::appender_output <io::formats::two_d::virtual_format> output_stream (file_name, 1, n, m, 1, 0, 0, 0, 0, 0);
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
				TSM_ASSERT_DELTA ("IO failure in virtual format output", io::virtual_dumps [file_name].index <double> ("test", i, j), init_copy [i * m + j], TEST_TINY);
			}
		}
	
		for (int i = 0; i < n * m; ++i) {
			init [i] = 0.0;
		}
		
		io::formatted_input <io::formats::two_d::virtual_format> input_stream (file_name, n, m, 1, 0, 0, 0, 0, 0);
		
		input_stream.append <double> ("test", &init [0]);
		input_stream.from_file ();
			
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				TSM_ASSERT_DELTA ("IO failure in virtual format input", init [i * m + j], init_copy [i * m + j], TEST_TINY);
			}
		}
	}
};
