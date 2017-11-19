#pragma once

#include "util.h"

#include "../cw_lib/CVec.h"
#include "../cw_lib/CCamera.h"

#include <ios>
#include <vector>

template<class T>
std::ostream& operator<< (ostream& os, const CVec2<T>& p) {
	return os << p.x << " " << p.y;
}

template<class T>
std::ostream& operator<< (ostream& os, const CVec3<T>& p) {
	return os << p.x << " " << p.y << " " << p.z;
}

template<class T>
std::ostream& operator<< (ostream& os, const CVec4<T>& p) {
	return os << p.x << " " << p.y << " " << p.z << " " << p.w;
}

template<class T>
std::ostream& dumpVector(ostream& os, const vector<T>& p) {
	os << p.size();
	for (const auto& x : p) {
		os << " " << x;
	}
	return os;
}

template<class T>
std::istream& operator>>(std::istream& is, CVec2<T>& p) {
	return is >> p.x >> p.y;
}

template<class T>
std::istream& operator>>(std::istream& is, CVec3<T>& p) {
	return is >> p.x >> p.y >> p.z;
}

template<class T>
std::istream& operator>>(std::istream& is, CVec4<T>& p) {
	return is >> p.x >> p.y >> p.z >> p.w;
}

inline
std::vector<CVec2f> project(CCamera cam, const std::vector<CVec3f>& points) {
	std::vector<CVec2f> result(points.size());
	for (size_t i = 0; i < points.size(); ++i)
		cam.WorldToImgCoords(points[i], result[i]);
	return result;
}

template<class To = cv::Vec2f, class From>
std::vector<To> project(CCamera cam, const std::vector<From>& points) {
	std::vector<To> result;
	for (size_t i = 0; i < points.size(); ++i) {
		CVec2f tmp;
		cam.WorldToImgCoords({ (float)points[i][0], (float)points[i][1], (float)points[i][2] }, tmp);
		result.emplace_back(tmp.x, tmp.y);
	}
	return result;
}
