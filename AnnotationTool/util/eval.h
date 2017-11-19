#pragma once

#include <opencv2/core.hpp>

namespace eval {

	enum SamplingType {
		SAMPLING_RANGE_UNIFORM,
		SAMPLING_DATA_UNIFORM,
		TOTAL
	};

	double auc(const cv::Mat& greyMap, const cv::Mat& gt, int tolerance, cv::InputArray mask);

	double bestf(const cv::Mat& greyMap, const cv::Mat& gt, int tolerance, cv::InputArray mask, cv::OutputArray prec_recall = cv::noArray());

}