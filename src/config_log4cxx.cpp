/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"
#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

log4cxx::LoggerPtr logger = log4cxx::Logger::getLogger ("log");
int severity = 2;

log_config log_config_instance;

log4cxx::LevelPtr int_to_severity (int severity_index) {
	switch (severity_index) {
		case 0:
		return log4cxx::Level::getTrace (); break;
		case 1:
		return log4cxx::Level::getDebug (); break;
		case 2:
		return log4cxx::Level::getInfo (); break;
		case 3:
		return log4cxx::Level::getWarn (); break;
		case 4:
		return log4cxx::Level::getError (); break;
		default:
		return log4cxx::Level::getFatal (); break;
	}
}

log_config::log_config () {
	logger->setLevel (int_to_severity (severity));
}

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
	
	std::ostringstream convert;
	convert << id;
	logger = log4cxx::Logger::getLogger ("process " + convert.str ());
	logger->setLevel (int_to_severity (severity));
	logger->addAppender (new log4cxx::ConsoleAppender (new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n")));

	logger->addAppender (new log4cxx::FileAppender (new log4cxx::PatternLayout ("%d %-5p %c{2}: %C (%M %L) - %m%n"), "process_" + convert.str () + ".log", false));
}

void log_config::trace (std::stringstream stream) {
	LOG4CXX_TRACE (logger, stream);
}

void log_config::debug (std::stringstream stream) {
	LOG4CXX_DEBUG (logger, stream);
}

void log_config::info (std::stringstream stream) {
	LOG4CXX_INFO (logger, stream);
}

void log_config::warn (std::stringstream stream) {
	LOG4CXX_WARN (logger, stream);
}

void log_config::error (std::stringstream stream) {
	LOG4CXX_ERROR (logger, stream);
}

void log_config::fatal (std::stringstream stream) {
	LOG4CXX_FATAL (logger, stream);
}

