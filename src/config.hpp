/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#ifndef _LOG4CXX

#include <iostream>

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) if(log_config::severity==0){std::cout<<"TRACE "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) if(log_config::severity<=1){std::cout<<"DEBUG "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) if(log_config::severity<=2){std::cout<<"INFO "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#define WARN(str) if(log_config::severity<=3){std::cout<<"WARN "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) if(log_config::severity<=4){std::cout<<"ERROR "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) if(log_config::severity<=5){std::cout<<"FATAL "<<str<<std::endl;}else{}

#else

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) LOG4CXX_TRACE(log_config::logger,str)
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) LOG4CXX_DEBUG(log_config::logger,str)
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) LOG4CXX_INFO(log_config::logger,str)
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#define WARN(str) LOG4CXX_WARN(log_config::logger,str)
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) LOG4CXX_ERROR(log_config::logger,str)
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) LOG4CXX_FATAL(log_config::logger,str)

#endif // _LOG4CXX

#ifdef _LOG4CXX

#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

#endif // _LOG4CXX

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class log_config {
public:
	log_config () {
		
#ifdef _LOG4CXX
		
		logger->setLevel (int_to_severity (severity));
		logger->addAppender (new log4cxx::ConsoleAppender (new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n")));
		
#endif // _LOG4CXX
	}
	
	/*!**********************************************************************
	 * \brief Updates the severity of the logger
	 * 
	 * The severity scheme is as follows:
	 *  0: Trace: Indicates a code location, used for optimization
	 *  1: Debug: Used to diagnose bugs
	 *  2: Info: Produces relevant information to user (default)
	 *  3: Warn: May indicate an issue, depending on user intent
	 *  4: Error: A problem that can be handled has occurred
	 *  5: Fatal: A problem that cannot be handled has occurred 
	 * 
	 * \param severity_index The integer severity below which log statements are suppressed
	 ************************************************************************/
	static void update_severity (int severity_index);
	
	/*!**********************************************************************
	 * \brief Updates the name of the logger and log file
	 * 
	 * The name of the logger becomes element_id, and the file name 
	 * element_id.log.
	 * 
	 * \param id The integer id with which to label the log.
	 ************************************************************************/
	static void update_name (int id);

	static int severity; //!< The integer level of severity to be output in the log
	
#ifdef _LOG4CXX
	/*!*******************************************************************
	 * \brief Converts an integer severity to a severity LevelPtr
	 * 
	 * \return The corresponding LevelPtr object for the logger
	 *********************************************************************/
	static log4cxx::LevelPtr int_to_severity (int severity_index);
	
	static log4cxx::LoggerPtr logger; //!< A pointer to a log4cxx logger object
	
#endif // _LOG4CXX

};

#endif /* end of include guard: CONFIG_H_95CWOMPS */
