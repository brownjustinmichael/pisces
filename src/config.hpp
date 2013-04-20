// 
//! \file config.hpp
//  src
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#ifdef __APPLE__

#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
	
/*!*******************************************************************
 * \def TRACE(str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) LOG4CXX_TRACE(config::logger,str)
/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) LOG4CXX_DEBUG(config::logger,str)
/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) LOG4CXX_INFO(config::logger,str)
/*!*******************************************************************
 * \def WARN(str)
 * Logs a warn-level logging statement
 *********************************************************************/
#define WARN(str) LOG4CXX_WARN(config::logger,str)
/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) LOG4CXX_ERROR(config::logger,str)
/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) LOG4CXX_FATAL(config::logger,str)

#else // __APPLE__

/*!*******************************************************************
 * \def TRACE(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define TRACE(str)
/*!*******************************************************************
 * \def DEBUG(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define DEBUG(str)
/*!*******************************************************************
 * \def INFO(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define INFO(str)
/*!*******************************************************************
 * \def WARN(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define WARN(str)
/*!*******************************************************************
 * \def ERROR(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define ERROR(str)
/*!*******************************************************************
 * \def FATAL(str)
 * To be overwritten with logging implementation
 *********************************************************************/
#define FATAL(str)

#endif // __APPLE__

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class config {
public:
	static int severity; //!< The integer level of severity to be output in the log
	
#ifdef __APPLE__
	/*!*******************************************************************
	 * \brief Converts an integer severity to a severity LevelPtr
	 * 
	 * \return The corresponding LevelPtr object for the logger
	 *********************************************************************/
	static log4cxx::LevelPtr int_to_severity (int severity_index);
	
	static log4cxx::LoggerPtr logger; //!< The logger object
		
#endif // __APPLE__
};


#endif /* end of include guard: CONFIG_H_95CWOMPS */
