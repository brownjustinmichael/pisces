/*!***********************************************************************
 * \file logger.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "logger.hpp"

#ifdef _LOG4CPLUS
// If LOG4CPLUS is defined, SCons has found a usable instance of log4cplus, and that will be used for logging

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/layout.h>

namespace logger
{
	log4cplus::BasicConfigurator config; //!< The log4cplus configurator object, must be configured at startup
	log4cplus::Logger logger = log4cplus::Logger::getRoot (); //!< The log4cplus logger object from which the log statements are derived
	log4cplus::SharedAppenderPtr append; //!< A log4cplus appender object that handles the actual IO
	int severity = 2; //!< The severity of output

	log_config log_config_instance; //!< The single necessary config instance only needed to insure that config.configure () is called
	
	/*!**********************************************************************
	 * \brief Convert from the integer severity to a log4cplus LogLevel object
	 * 
	 * \return A log4cplus LogLevel associated with the given severity
	 ************************************************************************/
	log4cplus::LogLevel int_to_severity (int severity_index) {
		switch (severity_index) {
			case 0:
			return log4cplus::TRACE_LOG_LEVEL; break;
			case 1:
			return log4cplus::DEBUG_LOG_LEVEL; break;
			case 2:
			return log4cplus::INFO_LOG_LEVEL; break;
			case 3:
			return log4cplus::WARN_LOG_LEVEL; break;
			case 4:
			return log4cplus::ERROR_LOG_LEVEL; break;
			default:
			return log4cplus::FATAL_LOG_LEVEL; break;
		}
	}

	log_config::log_config () {
		// Configure the config object to enable loggin
		config.configure();
		
		// Set the log level to the default
		logger.setLogLevel (int_to_severity (severity));
	}
	
	void log_config::configure (int* argc, char*** argv, int id, std::string log_file) {
		// Read the command line arguments
		for (int i = 0; i < *argc; ++i) {
			// Look for the severity logging level, e.g. -D3
			if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'D')) {
				severity = atoi (&((*argv) [i] [2]));
				--*argc;
				for (int j = i; j < *argc; ++j) {
					(*argv) [j] = (*argv) [j + 1];
				}
				--i;
			}
			
			// Look for a log file name e.g. -L "logger.log"
			if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'L')) {
				log_file = (*argv) [i + 1];
				*argc -= 2;
				for (int j = i; j < *argc; ++j) {
					(*argv) [j] = (*argv) [j + 2];
					(*argv) [j + 1] = (*argv) [j + 3];
				}
				i -= 2;
			}
		}
		
		// Set up the log file
		std::ostringstream convert;
		convert << id;
	    logger.setLogLevel (int_to_severity (severity));
		if (log_file != "") {
			char buffer [log_file.size () * 2];
			snprintf (buffer, log_file.size () * 2, log_file.c_str (), id);
			append = new log4cplus::FileAppender (buffer);
			append->setLayout (std::auto_ptr<log4cplus::Layout> (new log4cplus::PatternLayout ("%d %-5p: (%M %L) - %m%n")));
			logger.addAppender (append);
		}
		
		// Set up logging to the terminal
		append->setLayout (std::auto_ptr<log4cplus::Layout> (new log4cplus::PatternLayout ("%d %-5p: (%M %L) - %m%n")));
		logger.addAppender (append);
	}
	
	void log_config::set_severity (int i_severity) {
		severity = i_severity;
	    logger.setLogLevel (int_to_severity (i_severity));
	}
} /* logger */

#else
// If LOG4CPLUS is not defined, just use stdout

namespace logger
{
	int severity = 2; //!< The severity of output

	log_config log_config_instance; //!< The single necessary config instance only needed to insure that config.configure () is called

	log_config::log_config () {}

	void log_config::configure (int* argc, char*** argv, int id, std::string log_file) {
		for (int i = 0; i < *argc; ++i) {
			if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'D')) {
				severity = atoi (&((*argv) [i] [2]));
				--*argc;
				for (int j = i; j < *argc; ++j) {
					(*argv) [j] = (*argv) [j + 1];
				}
			}
		}
	}
	
	void log_config::set_severity (i_severity) {
		severity = i_severity;
	}
} /* logger */

#endif
