/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "logger.hpp"

#ifdef _LOG4CPLUS

#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/layout.h>

namespace logger
{
	log4cplus::BasicConfigurator config;
	log4cplus::Logger logger = log4cplus::Logger::getRoot ();
	log4cplus::SharedAppenderPtr append;
	int severity = 2;

	log_config log_config_instance;

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
		config.configure();
	    logger.setLogLevel (int_to_severity (severity));
	}

	void log_config::configure (int* argc, char*** argv, int id, std::string log_file) {
		for (int i = 0; i < *argc; ++i) {
			if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'D')) {
				severity = atoi (&((*argv) [i] [2]));
				--*argc;
				for (int j = i; j < *argc; ++j) {
					(*argv) [j] = (*argv) [j + 1];
				}
				--i;
			}
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
		append->setLayout (std::auto_ptr<log4cplus::Layout> (new log4cplus::PatternLayout ("%d %-5p: (%M %L) - %m%n")));
		logger.addAppender (append);
	}
} /* logger */

#else

namespace logger
{
	int severity = 2;
	log_config log_config_instance;

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
} /* logger */

#endif
