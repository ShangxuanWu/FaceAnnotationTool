// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu

#include "debug_tool.h"

int main(int argc, char* argv[]) {
	// ********************************************************************
	// test print_mat
	cv::Mat matrix_float(3, 3, CV_32F);
	//std::cout << (matrix.type() == CV_32FC2 || CV_32FC1) << std::endl;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			matrix_float.at<float>(i, j) = 1;
		}
	}
	print_mat(matrix_float);

	cv::Mat matrix_double(3, 3, CV_64F);
	//std::cout << (matrix.type() == CV_32FC2 || CV_32FC1) << std::endl;

	for (int i = 0; i < matrix_double.rows; i++) {
		for (int j = 0; j < matrix_double.cols; j++) {
			matrix_double.at<double>(i, j) = 1;
		}
	}
	print_mat(matrix_double);


	// ********************************************************************
	// test print_vec_pts3d
	std::vector<cv::Point3d> pts_double;
	pts_double.push_back(cv::Point3d(1, 2, 3));
	pts_double.push_back(cv::Point3d(4, 5, 6));
	print_vec_pts3d(pts_double);

	std::vector<cv::Point3f> pts_float;
	pts_float.push_back(cv::Point3f(1, 2, 3));
	pts_float.push_back(cv::Point3f(4, 5, 6));
	print_vec_pts3d(pts_float);

	std::vector<cv::Point3i> pts_int;
	pts_int.push_back(cv::Point3i(1, 2, 3));
	pts_int.push_back(cv::Point3i(4, 5, 6));
	print_vec_pts3d(pts_int);


	// ********************************************************************
	// test print_vec_pts2d
	std::vector<cv::Point2d> pts_double_2d;
	pts_double_2d.push_back(cv::Point2d(1, 2));
	pts_double_2d.push_back(cv::Point2d(4, 5));
	print_vec_pts2d(pts_double_2d);

	std::vector<cv::Point2f> pts_float_2d;
	pts_float_2d.push_back(cv::Point2f(1, 2));
	pts_float_2d.push_back(cv::Point2f(4, 5));
	print_vec_pts2d(pts_float_2d);

	std::vector<cv::Point2i> pts_int_2d;
	pts_int_2d.push_back(cv::Point2i(1, 2));
	pts_int_2d.push_back(cv::Point2i(4, 5));
	print_vec_pts2d(pts_int_2d);


	// ********************************************************************
	// test print_vec
	std::vector<double> vec_double;
	vec_double.push_back(1);
	vec_double.push_back(2);
	print_vec(vec_double);

	std::vector<float> vec_float;
	vec_float.push_back(1);
	vec_float.push_back(2);
	print_vec(vec_float);

	std::vector<int> vec_int;
	vec_int.push_back(1);
	vec_int.push_back(2);
	print_vec(vec_int);


	std::cout << std::endl;
	std::cout << "Testing done! Everything is good............." << std::endl;
	system("pause");

	return 0;
}