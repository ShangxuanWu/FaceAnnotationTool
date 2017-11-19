#pragma once

#include "convolution3.h"

#include <opencv2/core.hpp>
#include <algorithm>

struct LineIterator3 {
	LineIterator3(const cv::Vec3i& from, const cv::Vec3i& to)
		: i(0)
		, from(from)
		, diff(to - from)
		, ticks(std::max(std::max(abs(diff(0)), abs(diff(1))), abs(diff(2))))
	{
	}

	cv::Vec3i operator*() {
		return from + diff * i / ticks;
	}

	void operator++() {
		++i;
	}

	operator bool() {
		return i <= ticks;
	}

	bool operator !() {
		return !(bool(*this));
	}

	int i;
	cv::Vec3i from, diff;
	int ticks;
};

template<class T>
void drawLine3(cv::Mat& done, const Vec3i& from, const Vec3i& to, const T& color, int rad) {
	int r2 = sqr(rad);
	for (LineIterator3 it(from, to); it; ++it) {
		Vec3i p = *it;
		for (int z = -rad; z <= rad; ++z) if (size_t(z + p[0]) < done.size[0]) {
			for (int y = -rad; y <= rad; ++y) if (sqr(z) + sqr(y) <= r2 && size_t(y + p[1]) < done.size[1]) {
				for (int x = -rad; x <= rad; ++x) if (sqr(z) + sqr(y) + sqr(x) <= r2 && size_t(x + p[2]) < done.size[2]) {
					done.at<T>(z + p[0], y + p[1], x + p[2]) = color;
				}
			}
		}
	}
}

void nonMax3(const cv::Mat& conf, const cv::Mat& dir, double sigma, cv::OutputArray nMax);