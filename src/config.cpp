/*!***********************************************************************
 * \file config.cpp
 * Spectral Element
 * 
 * Created by Justin Brown on 2013-04-23.
 * Copyright 2013 Justin Brown. All rights reserved.
 ************************************************************************/

#include "config.hpp"

#ifdef __APPLE__

log4cxx::LevelPtr config::int_to_severity (int severity_index) {
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

#endif // #ifdef __APPLE__