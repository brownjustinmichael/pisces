/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#include "config_log4cxx.hpp"
#include "config_log4cplus.hpp"
#include <iostream>
#include <string>

extern int severity;

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef TRACE
#define TRACE(str) if(severity<=0){std::cout<<"TRACE "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef DEBUG
#define DEBUG(str) if(severity<=1){std::cout<<"DEBUG "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#ifndef INFO
#define INFO(str) if(severity<=2){std::cout<<"INFO "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#ifndef WARN
#define WARN(str) if(severity<=3){std::cout<<"WARN "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#ifndef ERROR
#define ERROR(str) if(severity<=4){std::cout<<"ERROR "<<str<<std::endl;}else{}
#endif
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#ifndef FATAL
#define FATAL(str) if(severity<=5){std::cout<<"FATAL "<<str<<std::endl;}else{}
#endif

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class log_config {
public:
	log_config ();
	
	virtual ~log_config () {
		// printf ("Destroying log config\n");
	}
	
	static void configure (int *argc, char ***argv, int id = 0, std::string log_file = "");
};

#endif /* end of include guard: CONFIG_H_95CWOMPS */
