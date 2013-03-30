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


class config {
public:
	static log4cxx::LoggerPtr logger;
	static log4cxx::LevelPtr levels [];
};

#endif /* end of include guard: CONFIG_H_95CWOMPS */
