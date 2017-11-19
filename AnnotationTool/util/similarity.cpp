#include "similarity.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>

#include <vector>

BOOST_GEOMETRY_REGISTER_POINT_3D(cv::Vec3d, double, boost::geometry::cs::cartesian, val[0], val[1], val[2])

BOOST_GEOMETRY_REGISTER_LINESTRING(std::vector<cv::Vec3d>)

double distanceAsym(const std::vector<cv::Vec3d>& s1, const std::vector<cv::Vec3d>& s2) {
	double bestD = 0.0;
	for (size_t i = 0; i < s2.size(); ++i) {
		bestD += boost::geometry::distance(s2[i], s1);
	}
	return bestD / s2.size();
}
