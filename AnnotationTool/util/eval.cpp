#include "eval.h"
#include "util.h"

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <vector>
#include <stdexcept>
#include <iterator>
#include <algorithm>

using namespace eval;
using namespace cv;

std::vector<std::pair<float, bool>> sortData(const cv::Mat& pred, const cv::Mat& gtMask, cv::InputArray mask, size_t& P, size_t& N) {
	std::vector<cv::Point> points;
	cv::findNonZero(mask, points);

	std::vector<std::pair<float, bool>> thresholds;
	P = 0;
	N = 0;
	for (const cv::Point& p : points) {
		thresholds.emplace_back(pred.at<float>(p), gtMask.at<uchar>(p));
		if (gtMask.at<uchar>(p))
			++P;
		else
			++N;
	}
	std::sort(thresholds.begin(), thresholds.end());
	return thresholds;
}

double eval::auc(const cv::Mat& greyMap, const cv::Mat& gt, int tolerance, cv::InputArray mask) {
	Mat pred, gtMask = gt > 0, gtMaskD;
	cv::dilate(gtMask, gtMaskD, Mat::ones(tolerance, tolerance, CV_8UC1));

	greyMap.convertTo(pred, CV_32F);

	size_t P = 0;
	size_t N = 0;

	auto thresholds = sortData(pred, gtMask, mask.getMat() & (gtMask | ~gtMaskD), P, N);
	thresholds.emplace_back(INFINITY, true);

	size_t curTP = P;
	size_t prevTP = P;
	size_t curFP = N;
	size_t prevFP = N;
	
	double res = 0;
	float prev = -INFINITY;
	for (auto pair : thresholds) {
		if (pair.first != prev) {
			res += double(prevTP + curTP) / P * double(prevFP - curFP) / N / 2;
			prevTP = curTP;
			prevFP = curFP;
			prev = pair.first;
		}
		if (pair.second)
			--curTP;
		else
			--curFP;
	}
	return res;
}

double eval::bestf(const cv::Mat& greyMap, const cv::Mat& gt, int tolerance, cv::InputArray mask, cv::OutputArray prec_recall) {
	Mat pred, gtMask = gt > 0, gtMaskD;
	cv::dilate(gtMask, gtMaskD, Mat::ones(tolerance, tolerance, CV_8UC1));

	greyMap.convertTo(pred, CV_32F);

	size_t P = 0;
	size_t N = 0;

	auto thresholds = sortData(pred, gtMask, mask.getMat() & (gtMask | ~gtMaskD), P, N);
	thresholds.emplace_back(INFINITY, true);

	size_t curTP = P;
	size_t curFP = N;

	bool collect = prec_recall.needed();
	std::vector<cv::Vec2f> pr;

	double bestF = 0;
	float prev = -INFINITY;
	for (auto pair : thresholds) {
		if (pair.first != prev) {
			double precision = double(curTP) / (curTP + curFP);
			double recall =    double(curTP) / P;
			double f = 2 * precision * recall / (precision + recall);
			relax<std::less>(bestF, f);
			if (collect)
				pr.emplace_back(precision, recall);
			prev = pair.first;
		}
		if (pair.second)
			--curTP;
		else
			--curFP;
	}

	if (collect) {
		prec_recall.create(pr.size(), 1, CV_32FC2);
		std::copy(pr.begin(), pr.end(), prec_recall.getMat().ptr<Vec2f>(0));
	}

	return bestF;
}
