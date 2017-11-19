#pragma once

#include <opencv2/core.hpp>

void convolve3(const cv::Mat& src, const cv::Mat& kernel, cv::OutputArray out);

void gaussianBlur3(const cv::Mat& src, double sigma, size_t kernelSize, cv::OutputArray out);

void dilate3(const cv::Mat& src, const cv::Mat& kernel, cv::OutputArray out);
