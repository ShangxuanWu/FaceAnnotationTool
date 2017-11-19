#include "util.h"
#include "cutil.h"

#include <Eigen/Core>

#include <fstream>

using namespace cv;

const char DIVIDERS[] = {
#if WIN32
	'/', '\\'
#else
	'/'
#endif
};

void report_range(Mat src) {
	Mat ch[3];
	split(src, ch);
	double a, b;
	for (int i = 0; i < src.channels(); ++i) {
		minMaxIdx(ch[i], &a, &b);
		printf("[%lf %lf] ", a, b);
	}
	printf("\n");
}

Mat to01(Mat a) {
	double min, max;
	minMaxLoc(a, &min, &max);
	return (a - min) / (max - min);
}

Mat angle2rgb(Mat angle, double scale, InputArray value) {
	Mat res;
	angle.convertTo(res, CV_8U, scale);
	Mat ch[3] = { res, res * 0 + 255, res * 0 + 255 };

	if (!value.empty()) {
		ch[2] = to01(value.getMat());
		ch[2].convertTo(ch[2], CV_8U, 255);
	}

	merge(ch, 3, res);
	cvtColor(res, res, CV_HSV2BGR);
	return res;
}

std::string replaceExt(const std::string &path, const std::string& ext) {
	size_t pos = path.find_last_of('.');
	return path.substr(0, pos) + ext;
}

std::string baseName(std::string fname) {
	for (char div : DIVIDERS) {
		size_t pos = fname.find_last_of(div);
		if (pos < fname.size()) {
			fname.erase(0, pos + 1);
		}
	}
	return fname;
}

std::string baseDir(std::string fname) {
	size_t rem = 0;
	for (char div : DIVIDERS)
	{
		size_t pos = fname.find_last_of(div);
		if (pos < fname.size()) {
			rem = max(rem, pos);
		}
	}
	if (rem)
		fname.erase(rem);
	else
		fname = "."; // no slashes in the path
	return fname;	
}

std::vector<std::vector<cv::Point>> read2DLines(const std::string& path) {
	std::vector<std::vector<cv::Point>> result;
	int n;
	for (std::ifstream ifs(path); ifs >> n;) {
		std::vector<cv::Point> tmp;
		for (int i = 0; i < n; ++i) {
			cv::Vec2f  p;
			ifs >> p[0] >> p[1];
			tmp.push_back(cv::Vec2i(p));
		}
		result.push_back(tmp);
	}
	return result;
}

cv::Mat mask2hmap(const cv::Mat& mask) {
	cv::Mat hmap;
	distanceTransform(mask == 0, hmap, CV_DIST_L2, CV_DIST_MASK_PRECISE);
	hmap = 1 / (1 + hmap);
	return hmap;
}

cv::Mat lines2mask(const std::vector<std::vector<cv::Point>>& lines, cv::Size size) {
	cv::Mat mask(size, CV_8UC1);
	mask = 0.;
	for (auto& pline : lines) {
		polylines(mask, pline, false, 255);
	}
	if (lines.empty()) {
		return cv::Mat();
	}
	return mask;
}

namespace details {
	template<int N, class T>
	std::ostream& dumpTrack(std::ostream& ofs, const T& track) {
		ofs << track.size();
		for (auto& p : track) {
			for (int j = 0; j < N; ++j)
				ofs << " " << p[j];
		}
		return ofs;
	}
}

#define DEFINE_DUMP_TRACK(T,N) template<> std::ostream& dumpTrack<T>(std::ostream& ofs, const std::vector<T>& track) { return details::dumpTrack<N>(ofs, track); }

#define DEFINE_DUMP_TRACK2(N) \
	DEFINE_DUMP_TRACK(cv::Vec ## N ## f, N) \
	DEFINE_DUMP_TRACK(CVec ## N ## f, N) \
	DEFINE_DUMP_TRACK(Eigen::Vector ## N ## f, N) \
	DEFINE_DUMP_TRACK(cv::Vec ## N ## d, N) \
	DEFINE_DUMP_TRACK(CVec ## N ## d, N) \
	DEFINE_DUMP_TRACK(Eigen::Vector ## N ## d, N) \

DEFINE_DUMP_TRACK2(2)
DEFINE_DUMP_TRACK2(3)

template<> std::ostream& dumpTrack<cv::Point>(std::ostream& ofs, const std::vector<cv::Point>& track) { 
	std::vector<cv::Vec<int, 2>> tmp(track.begin(), track.end());
	return details::dumpTrack<2>(ofs, tmp);
}
