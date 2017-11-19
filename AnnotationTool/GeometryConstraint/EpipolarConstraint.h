#ifndef EPIPOLAR_CONSTRAINT_H
#define EPIPOLAR_CONSTRAINT_H

#include <iostream>
#include <string>
#include <Windows.h>
#include <vector>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "EpipolarConstraint.h"
#include "../camera_propagation/pts_2d_conf.h"
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/camera_geometry.h"
#include "../camera_propagation/IdDataParser.h"

//void onMouse(int event, int x, int y);

void print_cv_mat(cv::Mat X);
__declspec(dllexport) cv::Mat getEpipolarLine(mycamera camera1, mycamera camera2, cv::Point original_point1);
__declspec(dllexport) cv::Mat getNullSpace(cv::Mat P);

class EpipolarConstraint {
public:
	EpipolarConstraint(std::string root_folder_);
	void getClickedPoint();
	void calculateEpipolarLine();
	void showEpipolarLine();
	void onMouse(int event, int x, int y);
	static void onMouse(int event, int x, int y, int, void* userdata);
private:
	void readCameraMatrix(std::string calibration_file_path, std::vector<mycamera>& camera_names);
	void calculateEpipolarLine(cv::Point original_point1);
	void getSubDirs(std::vector<std::string>& output, const std::string& path);
	//mycamera findCameraByName(std::string camera_name);
	int findCameraByName(std::string camera_name);

	std::string root_folder;
	std::vector<mycamera> cameras;
	std::vector<std::string> camera_names;
	int num1;
	int num2;
	int camera1_id;
	int camera2_id;
	float clicked_x;
	float clicked_y;

	cv::Mat image1;
	cv::Mat image2;
	cv::Mat image1_resized;
	cv::Mat image2_resized;

	int downsize_ratio = 8;
	int downsize_ratio2 = 8;
	cv::Size size;

	double line[3];
};


#endif // !EPIPOLAR_CONSTRAINT_H

