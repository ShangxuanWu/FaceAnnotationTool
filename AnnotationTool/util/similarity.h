#pragma once 

#include <opencv2/core.hpp>

#include <vector>

double distanceAsym(const std::vector<cv::Vec3d>& s1, const std::vector<cv::Vec3d>& s2);

inline double distanceSym(const std::vector<cv::Vec3d>& s1, const std::vector<cv::Vec3d>& s2) {
	return std::min(distanceAsym(s1, s2), distanceAsym(s1, s2));
}
