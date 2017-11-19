#include "convolution3.h"
#include "../util/util.h"
#include "../util/timer.h"

#include <opencv2/imgproc.hpp>

#include <iostream>

using namespace cv;

void gauss3(const cv::Mat& src, size_t ksize, double sigma, OutputArray out) {
	CV_Assert(src.dims == 3);
	CV_Assert(src.isContinuous());

	int dims[] = { src.size[0], src.size[1], src.size[2] };

	out.create(src.dims, dims, src.type());
	Mat dst = out.getMat();
	dst.setTo(0);

	Mat kernel(1, ksize, CV_64F);
	for (int i = 0; i < ksize; ++i) {
		kernel.at<double>(i) = (exp(-sqr((int)ksize/2 - i) / (2 * sqr(sigma))) / (sigma * sqrt(M_PI * 2)));
	}

	Mat tmp;
	for (int z = 0; z < dims[0]; ++z) {
		Mat subMat(2, dims + 1, src.type(), (void*)src.ptr(z, 0, 0));
		sepFilter2D(subMat, tmp, -1, kernel, kernel);
#pragma omp parallel for
		for (int i = max(0, z - (int)ksize / 2); i < min(z + (int)ksize / 2 + 1, dims[0]); ++i) {
			double weight = kernel.at<double>(ksize / 2 + z - i);
			Mat subDst(dst.size[1], dst.size[2], dst.type(), (void*)dst.ptr(i, 0, 0));
			subDst += tmp * weight;
		}
		if (z == 0) {
			double weight = 0;
			for (int j = ksize / 2; j > 0; --j) {
				weight += kernel.at<double>(ksize / 2 + j);
				if (j - 1 < dims[0]) {
					Mat subDst(dst.size[1], dst.size[2], dst.type(), (void*)dst.ptr(j - 1, 0, 0));
					subDst += tmp * weight;
				}
			}
		} 
		if (z == dims[0] - 1) {
			double weight = 0;
			for (int j = ksize / 2; j > 0; --j) {
				weight += kernel.at<double>(ksize / 2 + j);
				if (dims[0] >= j) {
					Mat subDst(dst.size[1], dst.size[2], dst.type(), (void*)dst.ptr(dims[0] - j, 0, 0));
					subDst += tmp * weight;
				}
			}
		}
	}
}

void convolve3(const cv::Mat& src, const cv::Mat& kernel, OutputArray out) {
	CV_Assert(src.dims == 3 && kernel.dims == 3);
	CV_Assert(src.isContinuous() && kernel.isContinuous());

	int dims[] = { src.size[0], src.size[1], src.size[2] };

	out.create(src.dims, dims, src.type());
	Mat dst = out.getMat();
	dst.setTo(0);

	Mat tmp;
	for (size_t i = 0; i < kernel.size[0]; ++i) {
		Mat subKernel(kernel.size[1], kernel.size[2], kernel.type(), (void*)kernel.ptr(i, 0, 0));
		for (int j = -(int)kernel.size[0] / 2; j < src.size[0] + kernel.size[0] / 2; ++j) {
			if ((size_t)j - i + kernel.size[0] / 2 < dst.size[0]) {
				int z = max(0, min((int)src.size[0] - 1, j));
				Mat subMat(src.size[1], src.size[2], src.type(), (void*)src.ptr(z, 0, 0));
				filter2D(subMat, tmp, -1, subKernel);
				Mat subDst(dst.size[1], dst.size[2], dst.type(), (void*)dst.ptr(j - i + kernel.size[0] / 2, 0, 0));
				subDst += tmp;
			}
		}
	}
}

void gaussianBlur3viaConvol(const cv::Mat& src, double sigma, size_t kernelSize, cv::OutputArray out) {
	double sigma2 = sigma*sigma;
	int dims[] = { kernelSize, kernelSize, kernelSize };
	cv::Mat kernel(3, dims, CV_64FC1);
	Vec3i anchor(kernelSize / 2, kernelSize / 2, kernelSize / 2);
	for (Vec3i pos(0, 0, 0); pos[2] < kernelSize; ++pos[2]) {
		for (pos[1] = 0; pos[1] < kernelSize; ++pos[1]) {
			for (pos[0] = 0; pos[0] < kernelSize; ++pos[0]) {
				Vec3d r = pos - anchor;
				kernel.at<double>(pos) = exp(-r.dot(r)/2/sigma2);
			}
		}
	}
	kernel /= pow(2 * M_PI * sigma2, 1.5);
	convolve3(src, kernel, out);
}

void gaussianBlur3(const cv::Mat& src, double sigma, size_t kernelSize, cv::OutputArray out) {
	gauss3(src, kernelSize, sigma, out);
}

void dilate3(const cv::Mat& src, const cv::Mat& kernel, OutputArray out) {
	CV_Assert(src.dims == 3 && kernel.dims == 3);
	CV_Assert(src.isContinuous() && kernel.isContinuous());

	int dims[] = { src.size[0], src.size[1], src.size[2] };

	out.create(src.dims, dims, src.type());
	Mat dst = out.getMat();
	dst.setTo(0);

	Mat tmp;
	for (size_t i = 0; i < kernel.size[0]; ++i) {
		Mat subKernel(kernel.size[1], kernel.size[2], kernel.type(), (void*)kernel.ptr(i, 0, 0));
		for (int j = 0; j < src.size[0]; ++j) {
			if ((size_t)j - i + kernel.size[0] / 2 < dst.size[0]) {
				Mat subMat(src.size[1], src.size[2], src.type(), (void*)src.ptr(j, 0, 0));
				dilate(subMat, tmp, subKernel);
				Mat subDst(dst.size[1], dst.size[2], dst.type(), (void*)dst.ptr(j - i + kernel.size[0] / 2, 0, 0));
				max(subDst, tmp, subDst);
			}
		}
	}
}
