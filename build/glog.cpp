#include <stdio.h>
#include <glog/logging.h>

int main(int argc, char* argv[]) {
	// Initialize Google's logging library.
	google::InitGoogleLogging(argv[0]);

	LOG (FATAL) << "Found " << " cookies";
}