#pragma once

#include "types.h"
#include "ceres_fwd.h"
#include "../cw_lib/CVec.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/segment.hpp>
#include <boost/geometry/geometries/register/linestring.hpp>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <opencv2/core.hpp>

#include <vector>

namespace bg = boost::geometry;

BOOST_GEOMETRY_REGISTER_POINT_3D(Eigen::Vector3d, double, bg::cs::cartesian, data()[0], data()[1], data()[2])
BOOST_GEOMETRY_REGISTER_POINT_3D(Eigen::Vector3f, float, bg::cs::cartesian, data()[0], data()[1], data()[2])
BOOST_GEOMETRY_REGISTER_POINT_3D(CVec3d, double, bg::cs::cartesian, x, y, z)
BOOST_GEOMETRY_REGISTER_POINT_3D(CVec3f, float, bg::cs::cartesian, x, y, z)
BOOST_GEOMETRY_REGISTER_POINT_3D(cv::Vec3d, double, bg::cs::cartesian, val[0], val[1], val[2])
BOOST_GEOMETRY_REGISTER_POINT_3D(cv::Vec3f, float, bg::cs::cartesian, val[0], val[1], val[2])

template<class T>
using Segment = std::pair<T, T>;

BOOST_GEOMETRY_REGISTER_SEGMENT_TEMPLATIZED(Segment, first, second)

BOOST_GEOMETRY_REGISTER_POINT_2D(Eigen::Vector2d, double, bg::cs::cartesian, data()[0], data()[1])
BOOST_GEOMETRY_REGISTER_POINT_2D(Eigen::Vector2f, float, bg::cs::cartesian, data()[0], data()[1])
BOOST_GEOMETRY_REGISTER_POINT_2D(CVec2d, double, bg::cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_POINT_2D(CVec2f, float, bg::cs::cartesian, x, y)
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Vec2d, double, bg::cs::cartesian, val[0], val[1])
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Vec2f, float, bg::cs::cartesian, val[0], val[1])
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Vec2i, int, bg::cs::cartesian, val[0], val[1])
BOOST_GEOMETRY_REGISTER_POINT_2D(cv::Point, int, bg::cs::cartesian, x, y)

BOOST_GEOMETRY_REGISTER_LINESTRING_TEMPLATED(std::vector)

template<class It1, class It2>
It2 subdivide(It1 begin, It1 end, double maxDist, It2 dst);

template<class It1, class It2>
It2 subdivide(It1 begin, It1 end, double maxDist, It2 dst) {
	auto prev = *begin++;
	auto last = prev;
	*dst++ = prev;
	bool lastAdded = true;
	for (; begin != end; ++begin) {
		auto& cur = *begin;
		auto diff = cur - prev;
		double dist = length(diff);
		int cnt = floor(dist / maxDist);
		for (int k = 1; k <= cnt; ++k)
			*dst++ = prev + diff * k / (cnt + 1);
		lastAdded = false;
		if (boost::geometry::distance(prev - cur) > 0.1 * maxDist) {
			*dst++ = cur;
			prev = cur;
			lastAdded = true;
		}
		last = cur;
	}
	if (!lastAdded) {
		*dst++ = last;
	}
	return dst;
}

namespace detail {

	template<class To>
	struct VectorCastBase {
		template<class... Args>
		static To create(const Args&... args) {
			typedef typename bg::coordinate_type<To>::type BaseType;
			using cv::saturate_cast;
			return To(saturate_cast<BaseType>(args)...);
		}
	};

	template<class To, class From, int D>
	struct VectorCast;

	template<class To, class From>
	struct VectorCast<To, From, 3> : VectorCastBase<To> {
		static To cast(const From& v) {
			return VectorCastBase<To>::create(bg::get<0>(v), bg::get<1>(v), bg::get<2>(v));
		}
	};

	template<class To, class From>
	struct VectorCast<To, From, 2> : VectorCastBase<To> {
		static To cast(const From& v) {
			return VectorCastBase<To>::create(bg::get<0>(v), bg::get<1>(v));
		}
	};

	template<class V>
	struct VectorCast<V, V, 2> {
		static const V& cast(const V& v) {
			return v;
		}
	};

	template<class V>
	struct VectorCast<V, V, 3> {
		static const V& cast(const V& v) {
			return v;
		}
	};
}

template<class To, class From>
typename std::enable_if<bg::dimension<To>::value == bg::dimension<From>::value, To>::type
vector_cast(const From& v) {
	return detail::VectorCast<To, From, bg::dimension<To>::value>::cast(v);
}

template<class Derived>
bool sameSide(const Eigen::MatrixBase<Derived>& p1, const Eigen::MatrixBase<Derived>& p2, const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<Derived>& b) {
	auto cp1 = (b - a).cross(p1 - a);
	auto cp2 = (b - a).cross(p2 - a);
	return cp1.dot(cp2) >= 0;
}

template<class Derived>
bool pointInTriangle(const Eigen::MatrixBase<Derived>& p, const Eigen::MatrixBase<Derived>& a, const Eigen::MatrixBase<Derived>& b, const Eigen::MatrixBase<Derived>& c) {
	return sameSide(p, a, b, c) && sameSide(p, b, a, c)
		&& sameSide(p, c, a, b);
}

template<class T>
typename T::Scalar pointTriangleDistance(const T& p, const T& v1, const T& v2, const T& v3) {
	typedef typename T::Scalar Scalar;
	T n = (v2 - v1).cross(v3 - v1);
	Eigen::Hyperplane<Scalar, T::SizeAtCompileTime> plane(n.normalized(), v1);
	T proj = plane.projection(p);
	if (pointInTriangle(proj, v1, v2, v3)) {
		return (proj - p).norm();
	}
	return std::min(std::min(
		bg::distance(p, Segment<T>(v1, v2)),
		bg::distance(p, Segment<T>(v1, v3))),
		bg::distance(p, Segment<T>(v2, v3)));
}

template<class T>
struct Rot {
	typedef Eigen::Matrix<T, 3, 1> Vec;

	Eigen::Matrix<T, 3, 3> rot;

	template<class V1, class V2>
	Rot(const V1& a, const V2& b)
		: rot(Eigen::Quaternion<T>::FromTwoVectors(vector_cast<Vec>(a), vector_cast<Vec>(b)))
	{
	}

	template<class V>
	V operator*(const V& x) const {
		return vector_cast<V>(Vec(rot * vector_cast<Vec>(x)));
	}

	template<class M>
	void matrix(M& matrix) const {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				matrix(i, j) = rot(i, j);
	}

	template<class M>
	M matrix() const {
		M m;
		matrix(m);
		return m;
	}
};

template<class T>
TrackT<T> normalize(const TrackT<T>& track, size_t count) {
	TrackT<T> newTrack;
	count = std::max<size_t>(2, count);
	double length = 0;
	for (size_t i = 0; i + 1 < track.size(); ++i)
		length += (asConst(track[i + 1]) - asConst(track[i])).norm();
	if (length < 1e-9)
		return newTrack;
	size_t j = 0;
	double before = 0, shiftLen = 0, segLen = (asConst(track[j + 1]) - asConst(track[j])).norm();
	for (size_t i = 0; i < count; ++i) {
		double next = length * i / (count - 1);
		while (before + shiftLen + 1e-12 < next) {
			if (before + segLen < next) {
				shiftLen = 0;
				before += segLen;
				++j;
				if (j == track.size() - 1)
					break;
				segLen = (asConst(track[j + 1]) - asConst(track[j])).norm();
			}
			else {
				shiftLen = next - before;
			}
		}
		newTrack.push_back(track[j]);
		if (shiftLen)
			newTrack.back() += (track[j + 1] - track[j]).normalized() * T(shiftLen);
	}
	return newTrack;
}

template<class T>
TracksT<T> normalize(TracksT<T> tracks, size_t count) {
	for (auto& track : tracks) {
		track = normalize(track, count);
	}
	return tracks;
}

template<class T, int N = 3>
void distance(const T point[N], const T q1[N], const T q2[N], T diff[N]) {
	typedef Eigen::Matrix<T, N, 1> Vec;
	typedef Eigen::Map<Vec> Map;
	typedef Eigen::Map<const Vec> CMap;

	Vec dq = CMap(q2) - CMap(q1);
	T t = dq.dot(CMap(point) - CMap(q2)) / dq.dot(dq);
	Map d(diff);
	if (t < 0.0) {
		d = Vec(CMap(point) - CMap(q1));
	}
	else if (t > 1.0) {
		d = CMap(point) - CMap(q1);
	}
	else {
		auto pPrime = CMap(q1) * t + CMap(q2) * (1.0 - t);
		d = CMap(point) - pPrime;
	}
}

template< int N = 3, class T, class V>
void distance(const T point[], const std::vector<V>& track, T dist[]) {
	CV_Assert(track.size() > 0);
	if (track.size() == 1) {
		Eigen::Map<Eigen::Matrix<T, N, 1>>(dist + 0) = track[0] - Eigen::Map<const Eigen::Matrix<T, N, 1>>(point);
		return;
	}
	double bestDistance = INFINITY;
	for (size_t i = 0; i + 1 < track.size(); ++i) {
		Eigen::Matrix<T, N, 1> diff;
		distance<T, N>(point, track[i].data(), track[i + 1].data(), diff.data());
		auto cDiff = diff.unaryExpr(std::ptr_fun(GetType<T>::removeJet));
		if (relax(bestDistance, cDiff.dot(cDiff)))
			Eigen::Map<Eigen::Matrix<T, N, 1>>(dist + 0) = diff;
	}
}
