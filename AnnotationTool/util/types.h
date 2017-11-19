#pragma once

#include <Eigen/Dense>

#include <string>
#include <vector>
#include <stdexcept>
#include <type_traits>

template<class T, int N = 3>
using VecT = Eigen::Matrix<T, N, 1>;

template<class T, int N = 3>
using MapT = typename std::conditional<std::is_const<T>::value,
	Eigen::Map<const VecT<typename std::remove_const<T>::type, N>>,
	Eigen::Map<VecT<T, N>>>::type;

template<class T, int N = 3>
using TrackT = std::vector<VecT<T, N>>;

template<class T, int N = 3>
using TracksT = std::vector<TrackT<T, N>>;

typedef VecT<double> Vec;
typedef MapT<double> Map;
typedef MapT<const double> CMap;
typedef TrackT<double> Track;
typedef TracksT<double> Tracks;
typedef std::vector<std::string> TracksData;

typedef VecT<double, 2> Vec2;
typedef MapT<double, 2> Map2;
typedef MapT<const double, 2> CMap2;
typedef TrackT<double, 2> Track2;
typedef TracksT<double, 2> Tracks2;

//template<class T>
//using CQMap = Eigen::Map<const Eigen::Quaternion<T>>;
//
template<class T>
using QMap = typename std::conditional<std::is_const<T>::value, 
	Eigen::Map<const Eigen::Quaternion<typename std::remove_const<T>::type>>, 
	Eigen::Map<Eigen::Quaternion<T>>>::type;

