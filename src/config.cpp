/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"

#ifdef _LOG4CXX

log4cxx::LevelPtr log_config::int_to_severity (int severity_index) {
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

#endif // _LOG4CXX

void log_config::update_severity (int severity_index) {
	severity = severity_index;

#ifdef _LOG4CXX
		logger->setLevel (int_to_severity (severity));
#endif // _LOG4CXX
}

void log_config::update_name (int id) {
#ifdef _LOG4CXX
		
	std::ostringstream convert;
	convert << id;
	logger = log4cxx::Logger::getLogger ("element_" + convert.str ());
	logger->setLevel (int_to_severity (severity));

	logger->addAppender (new log4cxx::ConsoleAppender (new log4cxx::PatternLayout ("%-5p %c{2}: %C (%M %L) - %m%n")));
	logger->addAppender (new log4cxx::FileAppender (new log4cxx::PatternLayout ("%d %-5p %c{2}: %C (%M %L) - %m%n"), "element_" + convert.str () + ".log", false));
		
#endif // _LOG4CXX
}

int log_config::severity = 2;

#ifdef _LOG4CXX

log4cxx::LoggerPtr log_config::logger = log4cxx::Logger::getLogger ("log");
log_config log_config_instance;

#endif // _LOG4CXX

