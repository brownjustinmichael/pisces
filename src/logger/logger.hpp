/*!***********************************************************************
 * \file logger.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#ifndef QUIET
/*!**********************************************************************
 * \def QUIET
 * During compile, logger will automatically suppress all computation of logging statements below the QUIET level. For production runs, QUIET should be at least 1, which will wipe all TRACE and DEBUG macros.
 ************************************************************************/
#define QUIET -1
#endif

#ifdef _LOG4CPLUS

#include <string>

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>

namespace logger
{
	extern log4cplus::Logger logger;
} /* logger */

/*!*******************************************************************
 * \def TRACE(str)
 * Logs a trace-level logging statement
 *********************************************************************/
#if QUIET < 0
#define TRACE(str) LOG4CPLUS_TRACE(logger::logger, str);
#else
#define TRACE(str)
#endif

/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement
 *********************************************************************/
#if QUIET < 1
#define DEBUG(str) LOG4CPLUS_DEBUG(logger::logger, str);
#else
#define DEBUG(str)
#endif

/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement
 *********************************************************************/
#if QUIET < 2
#define INFO(str) LOG4CPLUS_INFO(logger::logger, str);
#else
#define INFO(str)
#endif

/*!*******************************************************************
 * \def WARN(str)
 * Logs a warning-level logging statement
 *********************************************************************/
#if QUIET < 3
#define WARN(str) LOG4CPLUS_WARN(logger::logger, str);
#else
#define WARN(str)
#endif

/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement
 *********************************************************************/
#if QUIET < 4
#define ERROR(str) LOG4CPLUS_ERROR(logger::logger, str);
#else
#define ERROR(str)
#endif

/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#if QUIET < 5
#define FATAL(str) LOG4CPLUS_FATAL(logger::logger, str);
#else
#define FATAL(str)
#endif

#endif
#include <iostream>
#include <string>

namespace logger
{
	extern int severity;
} /* logger */


/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef TRACE
#if QUIET < 0
#define TRACE(str) if(logger::severity<=0){std::cout<<"TRACE - "<<str<<std::endl;}else{}
#else
#define TRACE(str)
#endif
#endif

/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef DEBUG
#if QUIET < 1
#define DEBUG(str) if(logger::severity<=1){std::cout<<"DEBUG - "<<str<<std::endl;}else{}
#else
#define DEBUG(str)
#endif
#endif

/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#ifndef INFO
#if QUIET < 2
#define INFO(str) if(logger::severity<=2){std::cout<<"INFO - "<<str<<std::endl;}else{}
#else
#define INFO(str)
#endif
#endif

/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#ifndef WARN
#if QUIET < 3
#define WARN(str) if(logger::severity<=3){std::cout<<"WARN - "<<str<<std::endl;}else{}
#else
#define WARN(str)
#endif
#endif

/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#ifndef ERROR
#if QUIET < 4
#define ERROR(str) if(logger::severity<=4){std::cout<<"ERROR - "<<str<<std::endl;}else{}
#else
#define ERROR(str)
#endif
#endif

/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#ifndef FATAL
#if QUIET < 5
#define FATAL(str) if(logger::severity<=5){std::cout<<"FATAL - "<<str<<std::endl;}else{}
#else
#define TRACE(str)
#endif
#endif

/*!*******************************************************************
 * \namespace logger
 * 
 * A namespace containing classes and objects relevant to log messages.
 *********************************************************************/
namespace logger
{
	/*!**********************************************************************
	 * \brief A config object associated with logging
	 * 
	 * All relevant contents of this object are static, so an instance need not be created by the user.
	 ************************************************************************/
	class log_config {
	public:
		log_config ();
		
		virtual ~log_config () {}
		
		/*!**********************************************************************
		 * \brief Configure the log object using information from the command line
		 * 
		 * \param argc The pointer to the number of command line arguments
		 * \param argv The pointer to the array of character arrays containing the command line arguments
		 * \param id The integer id for this logger (needed unique for parallel logging)
		 * \param log_file The string log file name, defaults to no log file
		 * 
		 * This function looks for a severity, given as an argument like -D3 for logging WARN and above, and for a log file name, given as -L "log_file.log". If no severity is given, it uses the default, 2: INFO and above. If no log file is given, it uses whatever is passed as the log_file argument. If log_file is "", nothing will be output to file.
		 ************************************************************************/
		static void configure (int *argc, char ***argv, int id = 0, std::string log_file = "");
		
		/*!**********************************************************************
		 * \brief Set the severity of the log messages
		 * 
		 * \brief severity The integer minimum severity to log (defaults to INFO)
		 * 
		 * The severities are as follows:
		 * 0 TRACE
		 * 1 DEBUG
		 * 2 INFO
		 * 3 WARN
		 * 4 ERROR
		 * 5 FATAL
		 * 
		 * For example, a severity of 2 will output all INFO, WARN, ERROR, and FATAL messages.
		 ************************************************************************/
		static void set_severity (int severity = 2);
	};
} /* logger */

#endif /* end of include guard: CONFIG_H_95CWOMPS */
