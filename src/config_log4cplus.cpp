/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/layout.h>

log4cplus::BasicConfigurator config;
log4cplus::Logger logger = log4cplus::Logger::getRoot ();
log4cplus::SharedAppenderPtr append;
int severity = 1;

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

log_config::log_config (int* argc, char*** argv, int id) {
	config.configure();
	
	while ((*argc > 1) && ((*argv) [1] [0] == '-')) {
		switch ((*argv) [1] [1]) {
			case 'D':
				severity = atoi (&((*argv) [1] [2]));
				break;
		}
		--*argc;
		++*argv;
	}
	std::ostringstream convert;
	convert << id;
    logger.setLogLevel (int_to_severity (severity));

	append = new log4cplus::FileAppender ("process_" + convert.str () + ".log");
	append->setLayout (std::auto_ptr<log4cplus::Layout> (new log4cplus::PatternLayout ("%d %-5p: (%M %L) - %m%n")));
	logger.addAppender (append);
}

