#pragma once

#include "eigen_io.h"
#include "util.h"
#include "types.h"
#include "../cw_lib/CVec.h"

#include <fstream>
#include <vector>
#include <string>

namespace txt_ndlines {
	template<class T>
	std::istream& read(std::istream& is, std::vector<std::vector<T>>& tracks, TracksData& data) {
		std::vector<T> track;
		std::string x;
		for (is;
			std::getline(readVector(is, track), x);
			tracks.push_back(std::move(track)), data.push_back(std::move(x)));
		return is;
	}

	template<class T>
	void read(const std::string& is, std::vector<std::vector<T>>& tracks, TracksData& data) {
		std::ifstream ifs(is);
		CV_Assert(ifs);
		read(ifs, tracks, data);
	}

	template<class T, class... Args>
	std::ostream& write(std::ostream& os, const std::vector<std::vector<T>>& tracks, const std::vector<Args>&... data) {
		for (size_t i = 0; i < tracks.size(); ++i) {
			dumpTrack(os, tracks[i], data[i]...) << std::endl;
		}
		return os;
	}

	template<class T, class... Args>
	void write(const std::string& os, const std::vector<std::vector<T>>& tracks, const std::vector<Args>&... data) {
		std::ofstream ofs(os);
		write(ofs, tracks, data...);
	}
}
