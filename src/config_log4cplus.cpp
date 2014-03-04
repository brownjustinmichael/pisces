/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include <iostream>
#include "config.hpp"
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <log4cplus/fileappender.h>
#include <log4cplus/layout.h>

log4cplus::BasicConfigurator config;
log4cplus::SharedAppenderPtr append;
log4cplus::Logger logger = log4cplus::Logger::getRoot ();

log_config log_config_instance;
int severity = 2;

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

void log_config::configure (int* argc, char*** argv, int id) {
	for (int i = 0; i < *argc; ++i) {
		if (((*argv) [i] [0] == '-') && ((*argv) [i] [1] == 'D')) {
			severity = atoi (&((*argv) [i] [2]));
			--*argc;
			for (int j = i; j < *argc; ++j) {
				(*argv) [j] = (*argv) [j + 1];
			}
		}
	}	

	std::ostringstream convert;
	convert << id;
    logger.setLogLevel (int_to_severity (severity));

	append = new log4cplus::FileAppender ("process_" + convert.str () + ".log");
	append->setLayout (std::auto_ptr<log4cplus::Layout> (new log4cplus::PatternLayout ("%d %-5p: (%M %L) - %m%n")));
	logger.addAppender (append);
}

void log_config::trace (std::stringstream &stream) {
	LOG4CPLUS_TRACE (logger, stream.str ());
}

void log_config::debug (std::stringstream &stream) {
	LOG4CPLUS_DEBUG (logger, stream.str ());
}

void log_config::info (std::stringstream &stream) {
	LOG4CPLUS_INFO (logger, stream.str ());
}

void log_config::warn (std::stringstream &stream) {
	LOG4CPLUS_WARN (logger, stream.str ());
}

void log_config::error (std::stringstream &stream) {
	LOG4CPLUS_ERROR (logger, stream.str ());
}

void log_config::fatal (std::stringstream &stream) {
	LOG4CPLUS_FATAL (logger, stream.str ());
}

