#pragma once

#include <opencv2/core.hpp>

void hessian2D(cv::Mat image, double sigma, cv::OutputArray h, bool level = false);

void frangiFilter2D(cv::Mat src, cv::OutputArray _dst, cv::OutputArray _ori, cv::OutputArray _scale, double sigmaStart, double sigmaStep, double sigmaEnd,
	double beta = 0.5, double c = 15.0, bool level = false, bool bw = true);
