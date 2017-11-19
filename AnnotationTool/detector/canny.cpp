#include "settings.h"
#include "../util/util.h"

using namespace cv;

void cannyThreshold(Mat conf, double lower, double upper, OutputArray res) {
	if (upper < lower)
		std::swap(upper, lower);
	Mat mask = conf >= upper;
	Mat notDone = conf >= lower;
	for (int i = 0; i < conf.rows; ++i) {
		for (int j = 0; j < conf.cols; ++j) {
			if (notDone.at<unsigned char>(i, j) == 255 && mask.at<char>(i, j)) {
				floodFill(notDone, Point(j, i), 127, NULL, 0, 0, 8);
			}
		}
	}
	Mat ans = notDone == 127;
	ans.copyTo(res);
}