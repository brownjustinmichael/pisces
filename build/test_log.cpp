#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <log4cplus/configurator.h>
#include <iomanip>
#include <iostream>
using namespace std;
using namespace log4cplus;
Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("main"));
void printMessages()
{
    LOG4CPLUS_TRACE(logger, "printMessages()");
    LOG4CPLUS_DEBUG(logger, "This is a DEBUG message");
    LOG4CPLUS_INFO(logger, "This is a INFO message");
    LOG4CPLUS_WARN(logger, "This is a WARN message");
    LOG4CPLUS_ERROR(logger, "This is a ERROR message");
    LOG4CPLUS_FATAL(logger, "This is a FATAL message");
}
int
main()
{
    initialize();
    BasicConfigurator config;
    config.configure();
    logger.setLogLevel(TRACE_LOG_LEVEL);
    cout << "*** calling printMessages() with TRACE set: ***" << endl;
    printMessages();
    logger.setLogLevel(DEBUG_LOG_LEVEL);
    cout << "\n*** calling printMessages() with DEBUG set: ***" << endl;
    printMessages();
    logger.setLogLevel(INFO_LOG_LEVEL);
    cout << "\n*** calling printMessages() with INFO set: ***" << endl;
    printMessages();
    logger.setLogLevel(WARN_LOG_LEVEL);
    cout << "\n*** calling printMessages() with WARN set: ***" << endl;
    printMessages();
    logger.setLogLevel(ERROR_LOG_LEVEL);
    cout << "\n*** calling printMessages() with ERROR set: ***" << endl;
    printMessages();
    logger.setLogLevel(FATAL_LOG_LEVEL);
    cout << "\n*** calling printMessages() with FATAL set: ***" << endl;
    printMessages();
    return 0;
}