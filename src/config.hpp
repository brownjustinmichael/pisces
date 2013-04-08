// 
//! \file config.hpp
//  src
//  
//  Created by Justin Brown on 2013-03-22.
//  Copyright 2013 Justin Brown. All rights reserved.
// 

#ifndef CONFIG_H_95CWOMPS
#define CONFIG_H_95CWOMPS

#include <log4cxx/logger.h>
#include <log4cxx/xml/domconfigurator.h>

#define TRACE(str) LOG4CXX_TRACE(config::logger,str)
#define DEBUG(str) LOG4CXX_DEBUG(config::logger,str)
#define INFO(str) LOG4CXX_INFO(config::logger,str)
#define WARN(str) LOG4CXX_WARN(config::logger,str)
#define ERROR(str) LOG4CXX_ERROR(config::logger,str)
#define FATAL(str) LOG4CXX_FATAL(config::logger,str)

class config {
public:
	static log4cxx::LoggerPtr logger;
	static log4cxx::LevelPtr levels [];
	static log4cxx::LevelPtr int_to_severity (int severity_index);
};

#endif /* end of include guard: CONFIG_H_95CWOMPS */
