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
 * Logs a trace-level logging statement
 *********************************************************************/
#define MTRACE(str) if(config::severity==0){std::cout<<"TRACE "<<str;}else{}
#define TRACE(int,str) if(config::severity==0&&int>=0){std::cout<<"element_"<<int<<" TRACE "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define MDEBUG(str) if(config::severity<=1){std::cout<<"DEBUG "<<str;}else{}
#define DEBUG(int,str) if(config::severity<=1&&int>=0){std::cout<<"element_"<<int<<" DEBUG "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement
 *********************************************************************/
#define MINFO(str) if(config::severity<=2){std::cout<<"INFO "<<str;}else{}
#define INFO(int,str) if(config::severity<=2&&int>=0){std::cout<<"element_"<<int<<" INFO "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def WARN(str)
 * Logs a warn-level logging statement
 *********************************************************************/
#define MWARN(str) if(config::severity<=3){std::cout<<"WARN "<<str;}else{}
#define WARN(int,str) if(config::severity<=3&&int>=0){std::cout<<"element_"<<int<<" WARN "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement
 *********************************************************************/
#define MERROR(str) if(config::severity<=4){std::cout<<"ERROR "<<str;}else{}
#define ERROR(int,str) if(config::severity<=4&&int>=0){std::cout<<"element_"<<int<<" ERROR "<<str<<std::endl;}else{}
/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define MFATAL(str) if(config::severity<=5){std::cout<<"FATAL "<<str;}else{}
#define FATAL(int,str) if(config::severity<=5&&int>=0){std::cout<<"element_"<<int<<" FATAL "<<str<<std::endl;}else{}

#else // __APPLE__

#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
	
/*!*******************************************************************
 * \def TRACE(str)
 * Logs a trace-level logging statement
 *********************************************************************/
#define MTRACE(str) LOG4CXX_TRACE(config::loggers[0],str)
#define TRACE(int,str) if(int>=0){LOG4CXX_TRACE(config::loggers[int],str)}else{}
/*!*******************************************************************
 * \def DEBUG(str)
 * Logs a debug-level logging statement
 *********************************************************************/
#define MDEBUG(str) LOG4CXX_DEBUG(config::loggers[0],str)
#define DEBUG(int,str) if(int>=0){LOG4CXX_DEBUG(config::loggers[int],str)}else{}
/*!*******************************************************************
 * \def INFO(str)
 * Logs an info-level logging statement
 *********************************************************************/
#define MINFO(str) LOG4CXX_INFO(config::loggers[0],str)
#define INFO(int,str) if(int>=0){LOG4CXX_INFO(config::loggers[int],str)}else{}
/*!*******************************************************************
 * \def WARN(str)
 * Logs a warn-level logging statement
 *********************************************************************/
#define MWARN(str) LOG4CXX_WARN(config::loggers[0],str)
#define WARN(int,str) if(int>=0){LOG4CXX_WARN(config::loggers[int],str)}else{}
/*!*******************************************************************
 * \def ERROR(str)
 * Logs an error-level logging statement
 *********************************************************************/
#define MERROR(str) LOG4CXX_ERROR(config::loggers[0],str)
#define ERROR(int,str) if(int>=0){LOG4CXX_ERROR(config::loggers[int],str)}else{}
/*!*******************************************************************
 * \def FATAL(str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#define MFATAL(str) LOG4CXX_FATAL(config::loggers[0],str)
#define FATAL(int,str) if(int>=0){LOG4CXX_FATAL(config::loggers[int],str)}else{}

#endif // __APPLE__

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class config {
public:
	static int severity; //!< The integer level of severity to be output in the log
	static int n_loggers;
	static int n_appenders;
	static std::string config_file;
	
#ifdef __APPLE__
	static std::vector<log4cxx::LoggerPtr> loggers;
	static std::vector<log4cxx::AppenderPtr> appenders;
	static log4cxx::LayoutPtr layout;
	static log4cxx::LayoutPtr console_layout;

	/*!*******************************************************************
	 * \brief Converts an integer severity to a severity LevelPtr
	 * 
	 * \return The corresponding LevelPtr object for the logger
	 *********************************************************************/
	static log4cxx::LevelPtr int_to_severity (int severity_index);
	
	static int make_main () {
		int i = n_loggers++;
		loggers.push_back (log4cxx::Logger::getLogger ("main"));
		loggers [i]->setLevel (config::int_to_severity (severity));
		
		int j = n_appenders++;
		appenders.push_back (new log4cxx::FileAppender (layout, "main.log", false));
		loggers [i]->addAppender (appenders [j]);
		
		j = n_appenders++;
		appenders.push_back (new log4cxx::ConsoleAppender (console_layout));
		loggers [i]->addAppender (appenders [j]);
		return i;
	}
	
	static int make_logger () {
		int i = n_loggers++;
		int j = n_appenders++;
		loggers.push_back (log4cxx::LoggerPtr (log4cxx::Logger::getLogger ("element_" + std::to_string (i))));
		loggers [i]->setLevel (config::int_to_severity (severity));
		
		appenders.push_back (new log4cxx::FileAppender (layout, "element_" + std::to_string (i) + ".log", false));
		loggers [i]->addAppender (appenders [j]);
		return i;
	}
	
#else
	
	static int make_main () {
		return 0;
	}
	
	static int make_logger () {
		++n_loggers;
		return n_loggers - 1; 
	}		

#endif // __APPLE__
};


#endif /* end of include guard: CONFIG_H_95CWOMPS */
