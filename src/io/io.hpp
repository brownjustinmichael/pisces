/*!***********************************************************************
 * \file io/io.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-08.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <vector>
#include <string>

#ifndef IO_HPP_C1E9B6EF
#define IO_HPP_C1E9B6EF

namespace io
{
	/*!*******************************************************************
	 * \brief A base class that is essentially a null header
	 * 
	 * Subclasses inheriting from this class should overload the 
	 * output_header method, which outputs a string that will represent the 
	 * header of the output file.
	 *********************************************************************/
	class header
	{
	public:
		virtual ~header () {}
		/*!*******************************************************************
		 * \brief Method that generates the header for an output object
		 * 
		 * This method should be overwritted by subclasses to be replaced with 
		 * more elaborate headers.
		 * 
		 * \return A string which will be used as the header in the associated output
		 *********************************************************************/
		virtual std::string output_header () {return std::string ("");}
	};
	
	/*!*******************************************************************
	 * \brief A simple header class
	 * 
	 * This class generates a simple header that contains just the names 
	 * of the columns.
	 *********************************************************************/
	class simple_header : public header
	{
	public:
		/*!*******************************************************************
		 * \param i_n_data_headers An integer number of column headers to be read from i_data_headers
		 * \param i_data_headers An array of strings that are the column headers
		 * \param i_comment_string The string to be used as the comment before the header lines
		 *********************************************************************/
		simple_header (int i_n_data_headers, std::string *i_data_headers, std::string i_comment_string = "#");
		
		virtual ~simple_header () {}
		
		/*!*******************************************************************
		 * \brief Make the header for the output file
		 * 
		 * \return The string to be streamed to the output file
		 *********************************************************************/
		std::string output_header () {return header_string;}

	private:
		std::string header_string; //!< The string form of the header
	};
	
	/*!*******************************************************************
	 * \brief An abstract output stream base class that generates output files
	 *********************************************************************/
	class output
	{
	public:
		/*!*******************************************************************
		 * \param i_header_ptr A pointer to the header object
		 * \param i_n The integer number of points in the data
		 * \param i_n_data_ptrs The integer number of pointers in i_data_ptrs
		 * \param i_data_ptrs A pointer to double pointers of data to output
		 *********************************************************************/
		output (header *i_header_ptr, int i_n, int i_n_data_ptrs = 0, double **i_data_ptrs = NULL);
		
		virtual ~output () {}
		
		/*!*******************************************************************
		 * \brief Append a double array to the list to be output
		 * 
		 * \param data_ptr A double pointer to the data to be the new column
		 *********************************************************************/
		virtual void append (double *data_ptr);
		
		/*!*******************************************************************
		 * \brief Append an integer array to the list to be output
		 * 
		 * \param data_ptr An integer pointer to the data to be the new column
		 *********************************************************************/
		virtual void append (int *data_ptr);
		
		/*!*******************************************************************
		 * \brief A virtual function to output the data to file
		 * 
		 * This function should be overwritten by subclasses, though it may 
		 * contain a call to this function, which will output with default 
		 * double representation in C++.
		 *********************************************************************/
		virtual void to_file () = 0;
		
	protected:
		int n; //!< An integer number of elements in each array
		int n_data_ptrs; //!< An integer number of arrays to output
		std::vector<double *> double_ptrs; //!< A vector of double pointers to the arrays of data (if an element is NULL, check int_ptrs at the same index instead)
		std::vector<int *> int_ptrs; //!< A vector of integer pointers to the arrays of data
		header *header_ptr; //!< A pointer to a header object, which contains the details regarding the construction of the header
		
		/*!*******************************************************************
		 * \brief Write all the tracked data to file using standard C conversions
		 * 
		 * \param file_name The string name of the file
		 *********************************************************************/
		void std_to_file (std::string file_name);
	};
	
	/*!*******************************************************************
	 * \brief A simple implementation of the output class
	 * 
	 * This class is a simple implementation of the output class.
	 *********************************************************************/
	class simple_output : public output
	{
	private:
		std::string file_name; //!< A string containing the file name where the class should output
	public:
		/*!*******************************************************************
		 * \param i_file_name A string of the file to which the data will be output
		 * \param i_n An integer number of data points contained within each array
		 * \param i_n_data_ptrs The number of double pointers in i_data_ptrs
		 * \param i_data_ptrs An array of double pointers to data arrays to be output
		 *********************************************************************/
		simple_output (std::string i_file_name, int i_n, int i_n_data_ptrs = 0, double **i_data_ptrs = NULL) : output (new header, i_n, i_n_data_ptrs, i_data_ptrs) {file_name = i_file_name;}
		
		/*!*******************************************************************
		 * \brief Outputs to file_name
		 *********************************************************************/
		void to_file () {std_to_file (file_name);} 
	};
	
	/*!*******************************************************************
	 * \brief An incremental implementation of the output class
	 * 
	 * This class is an implementation of the output class that increments 
	 * the file number with each iteration.
	 *********************************************************************/
	class incremental_output : public output
	{
	public:
		/*!*******************************************************************
		 * \param i_file_base A string containing the file name base (the string before the numbers, including the path)
		 * \param i_file_extension A string containing the file name extension (the string after the numbers, including the '.')
		 * \param i_int_width The total number of characters the integer in the file name can have
		 * \param i_header_ptr A pointer to the header object
		 * \param i_n The integer number of points in the data
		 * \param i_n_data_ptrs The number double pointers in i_data_ptrs
		 * \param i_data_ptrs An array of double pointers to data arrays to be output
		 *********************************************************************/
		incremental_output (std::string i_file_base, std::string i_file_extension, int i_int_width, header *i_header_ptr, int i_n, int i_n_data_ptrs = 0, double **i_data_ptrs = NULL) : output (i_header_ptr, i_n, i_n_data_ptrs, i_data_ptrs) {n_outputs = 0; int_width = i_int_width, file_base = i_file_base; file_extension = i_file_extension;}
		
		/*!*******************************************************************
		 * \brief Generates the next file name in the sequence
		 * 
		 * \return The incremented file name
		 *********************************************************************/
		std::string generate_file_name ();
		
		/*!*******************************************************************
		 * \brief Outputs to file
		 *********************************************************************/
		void to_file () {output::std_to_file (generate_file_name ());} 
		
	private:
		int n_outputs; //!< An incremental integer that will be appended to the end of each file_base
		int int_width; //!< The total number of characters the integer in the file name can have
		std::string file_base; //!< A string containing the file name base (the string before the numbers, including the path)
		std::string file_extension; //!< A string containing the file name extension (the string after the numbers, including the '.')
	};
} /* io */

#endif /* end of include guard: IO_HPP_C1E9B6EF */
