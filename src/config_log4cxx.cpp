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

log_config::log_config (int* argc, char*** argv, int id) {
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
	logger = log4cxx::Logger::getLogger ("process " + convert.str ());
	logger->setLevel (int_to_severity (severity));
	logger->addAppender (new log4cxx::ConsoleAppender (new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n")));

	logger->addAppender (new log4cxx::FileAppender (new log4cxx::PatternLayout ("%d %-5p %c{2}: %C (%M %L) - %m%n"), "process_" + convert.str () + ".log", false));
}

