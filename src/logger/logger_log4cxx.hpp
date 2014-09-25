/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_LOG4CXX_H_95CWOMPS
#define CONFIG_LOG4CXX_H_95CWOMPS

#ifdef _LOG4CXX

#include <log4cxx/logger.h>

extern log4cxx::LoggerPtr logger;

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define TRACE(str) LOG4CXX_TRACE(logger, str);
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define DEBUG(str) LOG4CXX_DEBUG(logger,str)
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#define INFO(str) LOG4CXX_INFO(logger,str)
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#define WARN(str) LOG4CXX_WARN(logger,str)
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#define ERROR(str) LOG4CXX_ERROR(logger,str)
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define FATAL(str) LOG4CXX_FATAL(logger,str)

#endif // _LOG4CXX

#endif /* end of include guard: CONFIG_H_95CWOMPS */
