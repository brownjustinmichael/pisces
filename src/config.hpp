/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#ifndef __APPLE__

#include <iostream>

/*!*******************************************************************
 * \def TRACE(str)
 * Logs a trace-level logging statement to the main log
 * \def TRACE(int,str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) if(log_config::severity==0){std::cout<<"TRACE "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement to the main log
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) if(log_config::severity<=1){std::cout<<"DEBUG "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement to the main log
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) if(log_config::severity<=2){std::cout<<"INFO "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def WARN(str)
 * Logs a warning-level logging statement to the main log
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#define WARN(str) if(log_config::severity<=3){std::cout<<"WARN "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement to the main log
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) if(log_config::severity<=4){std::cout<<"ERROR "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement to the main log
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) if(log_config::severity<=5){std::cout<<"FATAL "<<str<<std::endl;}else{}

#else

/*!*******************************************************************
 * \def TRACE(str)
 * Logs a trace-level logging statement to the main log
 * \def TRACE(int,str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) LOG4CXX_TRACE(log_config::logger,str)
/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement to the main log
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) LOG4CXX_DEBUG(log_config::logger,str)
/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement to the main log
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) LOG4CXX_INFO(log_config::logger,str)
/*!*******************************************************************
 * \def WARN(str)
 * Logs a warning-level logging statement to the main log
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#define WARN(str) LOG4CXX_WARN(log_config::logger,str)
/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement to the main log
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) LOG4CXX_ERROR(log_config::logger,str)
/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement to the main log
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) LOG4CXX_FATAL(log_config::logger,str)

#endif // __APPLE__

#ifdef __APPLE__

#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

#endif // __APPLE__

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class log_config {
public:
	log_config (int i_severity) {
		severity = i_severity;
		
#ifdef __APPLE__
		
		logger->setLevel (int_to_severity (severity));
		logger->addAppender (new log4cxx::ConsoleAppender (new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n")));
		
#endif // __APPLE__
	}
	
	~log_config () {}
	
	static int severity; //!< The integer level of severity to be output in the log
	
#ifdef __APPLE__
	static log4cxx::LoggerPtr logger;

	/*!*******************************************************************
	 * \brief Converts an integer severity to a severity LevelPtr
	 * 
	 * \return The corresponding LevelPtr object for the logger
	 *********************************************************************/
	static log4cxx::LevelPtr int_to_severity (int severity_index);

#endif // __APPLE__
	
	static void update_severity (int severity_index);
	
	static void update_name (int id);

};

#endif /* end of include guard: CONFIG_H_95CWOMPS */
