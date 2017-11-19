#pragma once

#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>

#include <fstream>

namespace cache {

	inline 
	bool newer(std::time_t refTime, const std::string& path) {
                if (refTime <= boost::filesystem::last_write_time(path))
                        return false;
                return true;
	}

	
	inline
	bool newer(std::time_t refTime, const boost::filesystem::path& path) {
		if (refTime <= boost::filesystem::last_write_time(path))
			return false;
		return true;
	}

	template<class Args>
        bool newer(std::time_t refTime, const std::vector<Args>& dependents) {
                for (const auto& path : dependents){
                        if (!newer(refTime, path))
                                return false;
                }
                return true;
        }
	

	template<class Args>
	size_t hash(const Args& dependents) {
		return boost::hash<Args>()(dependents);
	}

	template<class Args>
	std::ostream& writeHeader(std::ostream& os, const Args& dependents) {
		size_t h = hash(dependents);
		os.write((const char*)&h, sizeof(size_t));
		return os;
	}

	inline 
	std::istream& skipHeader(std::istream& is) {
		size_t h;
		is.read((char*)&h, sizeof(size_t));
		return is;
	}

	template<class T, class Args>
	bool needsUpdate(const T& cachePath, const Args& dependents) {
		if (boost::filesystem::exists(cachePath)) {
			size_t h;
			if (!std::ifstream(cachePath, std::ios::binary).read((char*)&h, sizeof(size_t)))
				return true;
			if (h != hash(dependents))
				return true;
			return !newer(boost::filesystem::last_write_time(cachePath), dependents);
		}
		return true;
	}

}
