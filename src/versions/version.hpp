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
#include <sstream>

/*!**********************************************************************
 * \namespace versions
 * 
 * \brief A namespace that contains the relevant information for any relevant versioning.
 ************************************************************************/
namespace versions
{
	/*!**********************************************************************
	 * \brief A struct that represents the version of a class
	 * 
	 * This class is designed to be used primarily through its inequality operator
	 ************************************************************************/
	struct version
	{
	private:
		int major; //!< The integer major version: large API changes and major physics additions
		int minor; //!< The integer minor version: small API changes, and large implementation changes
		int revision; //!< The integer revision version: small implementation changes and bug fixes
		int build; //!< The integer build: increment each time the code is built
		
	public:
		/*!**********************************************************************
		 * \brief Generate a version object from a string of the form "%d.%d.%d.%d"
		 * 
		 * \param versionStr A string of the form "%d.%d.%d.%d"
		 * 
		 * The string has the form major.minor.revision.build
		 ************************************************************************/
		version(std::string versionStr) {
			sscanf(versionStr.c_str(), "%d.%d.%d.%d", &major, &minor, &revision, &build);
		}
		
		/*!**********************************************************************
		 * \brief Compare two versions to see which is newer
		 ************************************************************************/
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
		
		/*!**********************************************************************
		 * \brief Compare two versions to see which is newer
		 ************************************************************************/
		bool operator<(const std::string &otherVersionString) const {
			version otherVersion (otherVersionString);
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
		
		/*!**********************************************************************
		 * \brief Compare two versions to see if they are the same
		 ************************************************************************/
		bool operator= (const version &otherVersion) const {
			return major == otherVersion.major && minor == otherVersion.minor && revision == otherVersion.revision && build == otherVersion.build;
		}

		/*!**********************************************************************
		 * \brief Get the string representation of the version
		 *
		 * \return The string representation of the version
		 ************************************************************************/
		operator std::string() {
			std::stringstream repr;
			repr << major << "." << minor << "." << revision << "." << build;
			return repr.str();
		}
	};
} /* versions */

#endif /* end of include guard: VERSION_HPP_9E6C97DA */
