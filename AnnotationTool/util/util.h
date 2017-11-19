#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
//#define _DEBUG_OUT

#include <opencv2/opencv.hpp>

#include <stddef.h>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <ios>
#include <iterator>
#include <functional>
#include <vector>

void report_range(cv::Mat src);

cv::Mat to01(cv::Mat a);

cv::Mat angle2rgb(cv::Mat angle, double scale = 180. / M_PI, cv::InputArray value = cv::noArray());

void nonMax(cv::Mat angle, cv::Mat src, double lambda, int directions, cv::OutputArray dst, int thickness = 2);

void cannyThreshold(cv::Mat conf, double lower, double upper, cv::OutputArray res);

// Writes float array into 16bit image (with loss)
void floatWrite(const std::string& destPath, const cv::Mat& data, cv::InputArray mask = cv::noArray());

cv::Mat floatRead(const std::string& path, int dstType);

cv::Mat unifiedRead(const std::string& path, int dstType);

cv::Mat mask2hmap(const cv::Mat& mat);

cv::Mat lines2mask(const std::vector<std::vector<cv::Point>>& path, cv::Size size);

std::vector<std::vector<cv::Point>> read2DLines(const std::string& path);

template<class T, int N>
bool inside(const cv::Mat& m, const cv::Vec<T, N>& pos);

/* If target > newVal then replace target with val and returns true 
   else returns false*/
template<class T, class Op = std::greater<T>>
bool relax(T& target, const T& newVal, Op op = std::greater<T>()) {
	if (op(target, newVal)) {
		target = newVal;
		return true;
	}
	return false;
}

/* If op(target, newVal) then replace target with val and returns true
 * else returns false. Here op is an object of template type Op<T>
 */
template<template<class> class Op, class T>
bool relax(T& target, const T& newVal, Op<T> op = Op<T>()) {
	return relax<T, Op<T>>(target, newVal, op);
}

template<class T>
auto sqr(const T& t) -> decltype(t * t) { return t * t; }

struct Longer {
	template<class T>
	bool operator()(const std::vector<T>& lhs, const std::vector<T>& rhs) const{
		return lhs.size() > rhs.size();
	}
};

#define REPORT_RANGE(x) do { printf(#x ": "); report_range(x); } while(false)

std::string replaceExt(const std::string &path, const std::string& ext);

std::string baseName(std::string fname);

std::string baseDir(std::string fname);

inline bool endsWith(const std::string& str, const std::string& suffix) {
	return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

template<class T>
struct Rollback {
	T& ref;
	const T val;

	Rollback(T& ref, const T& newVal)
		: ref(ref)
		, val(ref)
	{
		ref = newVal;
	}

	~Rollback() {
		ref = val;
	}
};

// implementation

template<int N> struct InsideImpl {
	template<class T, int M>
	static bool call(const int* p, const cv::Vec<T, M>& pos) {
		return 0 <= pos[N] && pos[N] < p[N] && InsideImpl<N - 1>::call(p, pos);
	}
};

template<> struct InsideImpl<-1> {
	template<class T, int M>
	static bool call(const int*, const cv::Vec<T, M>&) {
		return true;
	}
};

template<class T, int N>
bool inside(const int* p, const cv::Vec<T, N>& pos) {
	return InsideImpl<N - 1>::call(p, pos);
}

template<class T, int N>
bool inside(const cv::Mat::MSize& size, const cv::Vec<T, N>& pos) {
	return inside(size.p, pos);
}

template<class T, int N>
bool inside(const cv::Mat& m, const cv::Vec<T, N>& pos) {
	CV_Assert(m.dims >= N);
	return inside(m.size, pos);
}

inline cv::Mat unifiedRead(const std::string& path) {
	if (endsWith(path, ".yml"))
		return floatRead(path, CV_32F);
	return cv::imread(path, CV_LOAD_IMAGE_UNCHANGED);
}

template<class T>
std::istream& readVector(std::istream& is, std::vector<T>& p) {
	int sz;
	if (is >> sz) {
		p.resize(sz);
		for (auto& x : p) {
			is >> x;
		}
	}
	return is;
}

template<class T, int N>
std::istream& operator>>(std::istream& is, cv::Vec<T, N>& p) {
	for (int i = 0; i < N; ++i)
		is >> p[i];
	return is;
}

template<class T>
std::ostream& dumpTrack(std::ostream& ofs, const std::vector<T>& track);

template<class T, class... Args>
std::ostream& dumpTrack(std::ostream& ofs, const std::vector<T>& cur, Args&&... args) {
	dumpTrack(ofs, cur); 
	std::ostream* x[] = { &(ofs << " " << args)... };
	return ofs;
}
