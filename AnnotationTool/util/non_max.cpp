#include "../util/util.h"


using namespace cv;

void nonMaxHard(Mat angle, Mat src, double lambda, int directions, OutputArray dst) {
	Mat base = Mat::zeros(lambda * 2 + 1, lambda * 2 + 1, CV_8UC1);
	base.row(lambda).setTo(1);
	//double factor = tan(M_PI / directions * 2);
	//for (int i = 0; i < base.cols; ++i) {
	//	int y = abs(i - base.cols / 2) * factor;
	//	base.col(i).rowRange(base.rows / 2 - y, base.rows / 2 + y + 1).setTo(1);
	//}

	//base.colRange((base.cols - lambda) / 2, (base.cols + lambda) / 2).setTo(0);

#ifdef _DEBUG_OUT
	imwrite("base.png", base * 255);
#endif

	//#ifdef _DEBUG_OUT
	imwrite("conf.png", to01(src) * 255);
	//#endif

	Mat cur;
	for (int i = 0; i < directions; ++i) {

		double ang = i *  180. / directions;
		Mat rot = getRotationMatrix2D(Point(lambda, lambda), -ang, 1);
		Mat kernel;
		warpAffine(base, kernel, rot, base.size());
		kernel = kernel > 0;

#ifdef _DEBUG_OUT
		char buf[250];
		sprintf(buf, "k_%d.png", i);
		imwrite(buf, kernel);
#endif

		dilate(src, cur, kernel);

		src.setTo(0, angle == i & cur > src);
	}

	src.copyTo(dst);
}

#include <iostream>

void nonMaxBlur(Mat angle, Mat src, double lambda, int directions, OutputArray dst, int thickness) {

	if (lambda) {
		lambda /= 3;
		GaussianBlur(src, src, Size(int(lambda * 4) | 1, int(lambda * 4) | 1), lambda);
	}

#ifdef _DEBUG_OUT
	imwrite("blured.png", src);
#endif

	Point delta[] = { { 1, 0 }, { 1, 1 }, { 0, 1 }, { -1, 1 } };

	int dirMap[255];
	CV_Assert(directions < 255);
	CV_Assert(angle.type() == CV_32S);

	double range = (double)directions / 4;
	for (int i = 0; i < directions; ++i) {
		dirMap[i] = int(i / range + .5) % 4;
	}

	Mat mask(src.size(), CV_8UC1);
	mask.setTo(0);
	Rect rect(Point(), src.size());
	for (Point p; p.y < src.rows; ++p.y)
	for (p.x = 0; p.x < src.cols; ++p.x) {
		for (int dir = -1; dir <= 1; dir += 2) {
			Point neigh = p + thickness * dir * delta[dirMap[angle.at<int>(p)]];
			if (rect.contains(neigh) && src.at<float>(p) < src.at<float>(neigh)) {
				mask.at<char>(p) = 1;
				break;
			}
		}
	}
	src.setTo(0, mask);
#ifdef _DEBUG_OUT
	imwrite("non-max.png", src);
#endif
	src.copyTo(dst);
}

void nonMax(Mat angle, Mat src, double lambda, int directions, OutputArray dst, int thickness) {
	nonMaxBlur(angle, src, lambda, directions, dst, thickness);
}
