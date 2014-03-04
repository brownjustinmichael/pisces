/*!***********************************************************************
 * \file config.hpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#include <iostream>
#include <sstream>

/*!*******************************************************************
 * \brief A class containing the relevant configuration details, such as logging severity
 *********************************************************************/
class log_config {
public:
	log_config ();
	
	virtual ~log_config () {
		// printf ("Destroying log config\n");
	}
	
	void trace (std::stringstream &stream);

	void debug (std::stringstream &stream);

	void info (std::stringstream &stream);

	void warn (std::stringstream &stream);
	
	void error (std::stringstream &stream);
	
	void fatal (std::stringstream &stream);
	
	static void configure (int *argc, char ***argv, int id = 0);
};

extern log_config log_config_instance;
extern int severity;

/*!*******************************************************************
 * \def TRACE(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef TRACE
#define TRACE(str) if(severity<=0){std::stringstream stream; stream << str; log_config_instance.trace (stream);}else{}
#endif
/*!*******************************************************************
 * \def DEBUG(int,str)
 * Logs a debug-level logging statement
 *********************************************************************/
#ifndef DEBUG
#define DEBUG(str) if(severity<=1){std::stringstream stream; stream << str; log_config_instance.debug (stream);}else{}
#endif
/*!*******************************************************************
 * \def INFO(int,str)
 * Logs an info-level logging statement
 *********************************************************************/
#ifndef INFO
#define INFO(str) if(severity<=2){std::stringstream stream; stream << str; log_config_instance.info (stream);}else{}
#endif
/*!*******************************************************************
 * \def WARN(int,str)
 * Logs a warning-level logging statement
 *********************************************************************/
#ifndef WARN
#define WARN(str) if(severity<=3){std::stringstream stream; stream << str; log_config_instance.warn (stream);}else{}
#endif
/*!*******************************************************************
 * \def ERROR(int,str)
 * Logs an error-level logging statement
 *********************************************************************/
#ifndef ERROR
#define ERROR(str) if(severity<=4){std::stringstream stream; stream << str; log_config_instance.error (stream);}else{}
#endif
/*!*******************************************************************
 * \def FATAL(int,str)
 * Logs a fatal-level logging statement
 *********************************************************************/
#ifndef FATAL
#define FATAL(str) if(severity<=5){std::stringstream stream; stream << str; log_config_instance.fatal (stream);}else{}
#endif

#endif /* end of include guard: CONFIG_H_95CWOMPS */
