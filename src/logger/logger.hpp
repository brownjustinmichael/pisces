/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#ifndef QUIET
	#define QUIET -1
#endif

#include <string>

#ifdef _LOG4CPLUS

	#include <log4cplus/logger.h>
	#include <log4cplus/loggingmacros.h>

	namespace logger
	{
		extern log4cplus::Logger logger;
	} /* logger */

	/*!*******************************************************************
	 * \def TRACE(int,str)
	 * Logs a trace-level logging statement
	 *********************************************************************/
	#if QUIET >= 0
		#define TRACE(str)
	#else
		#define TRACE(str) LOG4CPLUS_TRACE(logger::logger, str);
	#endif
	/*!*******************************************************************
	 * \def DEBUG(int,str)
	 * Logs a debug-level logging statement
	 *********************************************************************/
	#if QUIET >= 1
		#define DEBUG(str)
	#else
		#define DEBUG(str) LOG4CPLUS_DEBUG(logger::logger, str);
	#endif
	/*!*******************************************************************
	 * \def INFO(int,str)
	 * Logs an info-level logging statement
	 *********************************************************************/
	#if QUIET >= 2
		#define INFO(str)
	#else
		#define INFO(str) LOG4CPLUS_INFO(logger::logger, str);
	#endif
	/*!*******************************************************************
	 * \def WARN(int,str)
	 * Logs a warning-level logging statement
	 *********************************************************************/
	#if QUIET >= 3
		#define WARN(str)
	#else
		#define WARN(str) LOG4CPLUS_WARN(logger::logger, str);
	#endif
	/*!*******************************************************************
	 * \def ERROR(int,str)
	 * Logs an error-level logging statement
	 *********************************************************************/
	#if QUIET >= 4
		#define ERROR(str)
	#else
		#define ERROR(str) LOG4CPLUS_ERROR(logger::logger, str);
	#endif
	/*!*******************************************************************
	 * \def FATAL(int,str)
	 * Logs a fatal-level logging statement
	 *********************************************************************/
	#if QUIET >= 5
		#define FATAL(str)
	#else
		#define FATAL(str) LOG4CPLUS_FATAL(logger::logger, str);
	#endif

#endif
	
#include <iostream>

namespace logger
{
	extern int severity;
} /* logger */

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef TRACE
#define TRACE(str) if(logger::severity<=0){std::cout<<"TRACE - "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef DEBUG
#define DEBUG(str) if(logger::severity<=1){std::cout<<"DEBUG - "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#ifndef INFO
#define INFO(str) if(logger::severity<=2){std::cout<<"INFO - "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#ifndef WARN
#define WARN(str) if(logger::severity<=3){std::cout<<"WARN - "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#ifndef ERROR
#define ERROR(str) if(logger::severity<=4){std::cout<<"ERROR - "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#ifndef FATAL
#define FATAL(str) if(logger::severity<=5){std::cout<<"FATAL - "<<str<<std::endl;}else{}
#endif

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
namespace logger
{
	class log_config {
	public:
		log_config ();

		virtual ~log_config () {
			// printf ("Destroying log config\n");
		}

		static void configure (int *argc, char ***argv, int id = 0, std::string log_file = "");

		static void set_severity (int severity = 2);
	};
} /* logger */

#endif /* end of include guard: CONFIG_H_95CWOMPS */
