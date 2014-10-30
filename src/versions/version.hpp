/*!**********************************************************************
 * \file version.hpp
 * /Users/justinbrown/Dropbox/pisces/build
 * 
 * Created by Justin Brown on 2014-10-30.
 * Copyright 2014 Justin Brown. All rights reserved.
 ************************************************************************/

#ifndef VERSION_HPP_9E6C97DA
#define VERSION_HPP_9E6C97DA

#include <string>

namespace versions
{
	struct version
	{
		version(std::string versionStr) {
			sscanf(versionStr.c_str(), "%d.%d.%d.%d", &major, &minor, &revision, &build);
		}

		bool operator<(const version &otherVersion) const {
			if (major < otherVersion.major)
				return true;
			if (major > otherVersion.major)
				return false;
			if (minor < otherVersion.minor)
				return true;
			if (minor > otherVersion.minor)
				return false;
			if (revision < otherVersion.revision)
				return true;
			if (revision > otherVersion.revision)
				return false;
			if (build < otherVersion.build)
				return true;
			return false;
		}
		
		bool operator= (const version &otherVersion) const {
			return major == otherVersion.major && minor == otherVersion.minor && revision == otherVersion.revision && build == otherVersion.build;
		}

		int major, minor, revision, build;
	};
} /* versions */

#endif /* end of include guard: VERSION_HPP_9E6C97DA */
