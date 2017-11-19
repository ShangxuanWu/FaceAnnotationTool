#include "util.h"

#include <opencv2/opencv.hpp>

using namespace cv;

void floatWrite(const std::string& destPath, const cv::Mat& data, InputArray mask) {
	double min = INFINITY, max = -INFINITY;
	if (mask.empty() || countNonZero(mask) > 0)
		minMaxIdx(data, &min, &max, NULL, NULL, mask);

#if _DEBUG_OUT
	std::cerr << "[" << id << "] data range: " << min << " " << max << std::endl;
#endif

	std::string fname = baseName(destPath) + ".png";
	std::string dir = baseDir(destPath);
	std::string dataPath = dir.empty() ? fname : format("%s/%s", dir.c_str(), fname.c_str());

	FileStorage fs(destPath, FileStorage::WRITE);

	if (min < INFINITY && min <  max && max > -INFINITY) {
		fs << "Min" << min << "Max" << max << "DataPath" << fname;
		Mat encoded;
		data.convertTo(encoded, CV_16U, USHRT_MAX / (max - min), -USHRT_MAX * min / (max - min));
		imwrite(dataPath, encoded);
	}
	else {
		fs << "Min" << 0 << "Max" << 0 << "DataPath" << fname;
		Mat encoded(data.size(), CV_16UC1);
		encoded = 0.0;
		imwrite(dataPath, encoded);
		std::cerr << "Empty data range: " << dataPath << std::endl;
	}
}

cv::Mat floatRead(const std::string& path, int dstType) {
	CV_Assert(dstType == CV_64F || dstType == CV_32F);

	std::string fname;
	double min, max;

	{
		FileStorage fs(path, FileStorage::READ);
		fs["Min"] >> min;
		fs["Max"] >> max;
		fs["DataPath"] >> fname;
	}
	
	fname = format("%s/%s", baseDir(path).c_str(), fname.c_str());
	Mat data = imread(fname, CV_LOAD_IMAGE_UNCHANGED);
	if(data.empty()) {
		throw std::runtime_error(format("Failed to load file: %s", fname.c_str()));
	}
	data.convertTo(data, dstType, (max - min) / USHRT_MAX, min);

	return data;
}
