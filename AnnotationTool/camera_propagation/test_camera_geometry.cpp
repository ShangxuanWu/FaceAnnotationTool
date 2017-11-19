// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu


#include "myheader.h"

#include <opencv2/highgui/highgui.hpp>

#include "debug_tool.h"
#include "mycamera.h"
#include "camera_geometry.h"
#include "IdDataParser.h"
#include "io_point.h"
#include "pts_2d_conf.h"
#include "pts_3d_conf.h"
#include "type_conversion.h"
#include "math_functions.h"


#define consider_dist_test		true
#define consider_skew_test		false
#define resize_factor_test		1
#define number_pts_load_test	68
#define num_ransac_test			100
//#define default_frame_test		0
#define seed_test				std::time(NULL)

// TODO: two view should be exact equal
// TODO: epipolar line should be exact precise

int main(int argc, char** argv)
{
	// read calibration file
	std::vector<mycamera> camera_cluster;
	camera_cluster.clear();
	LoadIdCalibration("calibration.txt", camera_cluster, consider_skew_test, consider_dist_test);		// consider distortion

	//std::cout << "a" == "a" << std::endl;

	// test get_two_random_camera_id
	//std::cout << std::endl;
	//std::cout << "Testing get_two_random_camera_id............." << std::endl;
	//std::map<std::string, pts_2d_conf> pts_2d_given_test;
	//pts_2d_given_test["test1"] = pts_2d_conf(2, 3, 1);
	//pts_2d_given_test["test2"] = pts_2d_conf(2, 3, 0.1);
	//pts_2d_given_test["test3"] = pts_2d_conf(2, 3, 1);
	//pts_2d_given_test["test4"] = pts_2d_conf(2, 3, 0);
	//std::vector<int> view_selected_random(2);
	//std::srand(seed_test);
	//for (int i = 0; i < 20; i++) {
	//	get_two_random_camera_id(pts_2d_given_test, view_selected_random, std::rand());
	//	print_vec(view_selected_random);
	//}

	// test multi_view_projection
	//std::cout << std::endl;
	//std::cout << "Testing multi_view_projection............." << std::endl;
	//std::vector<cv::Point3d> pts_src_projection;
	//pts_src_projection.push_back(cv::Point3d(1, 2, 500));
	////pts_src_projection.push_back(cv::Point3d(2, 4, 500));
	//std::map<std::string, std::vector<cv::Point2d>> pts_dst_projection;
	//multi_view_projection(pts_src_projection, camera_cluster, pts_dst_projection, consider_dist);
	//print_vec_pts2d(pts_dst_projection["330006"]);

	//mycamera camera_projection;
	//for (int i = 0; i < camera_cluster.size(); i++) {
	//	if (strcmp(camera_cluster[i].name.c_str(), "330006") == 0)
	//		camera_projection = camera_cluster[i];
	//}
	//print_mat(camera_projection.getProjectionMatrix());
	//std::vector<double> pts_src_vec = {1, 2, 500, 1};
	//cv::Mat pts_src_mat = camera_projection.getProjectionMatrix() * cv::Mat(pts_src_vec);
	//print_mat_info(pts_src_mat);
	//std::vector<double> pts_dst_vec = mat2vec(pts_src_mat);
	//for (int i = 0; i < pts_dst_vec.size() - 1; i++)
	//	pts_dst_vec[i] = pts_dst_vec[i] / pts_dst_vec[pts_dst_vec.size() - 1];
	//print_vec(pts_dst_vec);


	// test triangulation_from_two_views and undistort_single_point
	std::cout << std::endl;
	std::cout << "Testing triangulation_from_two_views and undistort_single_point............." << std::endl;
	std::vector<pts_2d_conf> pts_src1, pts_src2;
	load_points_with_conf("test\\cam330006\\00000.pose", pts_src1, resize_factor_test, number_pts_load_test);
	load_points_with_conf("test\\cam330014\\00000.pose", pts_src2, resize_factor_test, number_pts_load_test);
	print_vec_pts_2d_conf(pts_src1);
	print_vec_pts_2d_conf(pts_src2);
	std::vector<pts_3d_conf> pts_3d_two_view_new;
	std::vector<cv::Point3d> pts_3d_two_view_old;
	mycamera camera1, camera2;
	for (int i = 0; i < camera_cluster.size(); i++) {
		if (strcmp(camera_cluster[i].name.c_str(), "330006") == 0)
			camera1 = camera_cluster[i];
		if (strcmp(camera_cluster[i].name.c_str(), "330014") == 0)
			camera2 = camera_cluster[i];
	}
	triangulation_from_two_views(pts_src1, pts_src2, camera1, camera2, pts_3d_two_view_new, consider_dist_test);
	triangulation_from_two_views(conf2cv_vec_pts2d(pts_src1), conf2cv_vec_pts2d(pts_src2), camera1, camera2, pts_3d_two_view_old, consider_dist_test);
	std::map<std::string, std::vector<cv::Point2d>> pts_dst_multiview_new, pts_dst_multiview_old;
	multi_view_projection(conf2cv_vec_pts3d(pts_3d_two_view_new), camera_cluster, pts_dst_multiview_new, consider_dist_test);
	multi_view_projection(pts_3d_two_view_old, camera_cluster, pts_dst_multiview_old, consider_dist_test);
	print_vec_pts3d(pts_3d_two_view_old);
	print_vec_pts2d(pts_dst_multiview_old["330006"]);
	print_vec_pts2d(pts_dst_multiview_old["330014"]);
	
	print_vec_pts_3d_conf(pts_3d_two_view_new);
	print_vec_pts2d(pts_dst_multiview_new["330006"]);
	print_vec_pts2d(pts_dst_multiview_new["330014"]);
	std::cout << "Please check if they are approximately equal." << std::endl;


	/* test multiview_optimization_multiple_points, multiview_optimization_single_point, multiview_ransac_single_point, calculate_projection_error */
	std::cout << std::endl;
	std::cout << "Testing multiview_optimization_multiple_point, multiview_optimization_single_point, multiview_ransac_multiple_points, \
	multiview_ransac_single_point, calculate_projection_error............." << std::endl;
	int frame = 0;
	std::map<std::string, std::vector<pts_2d_conf>> pts_src_optimization;
	load_points_with_conf_multiview("image_original_conf0", camera_cluster, frame, pts_src_optimization, resize_factor_test, number_pts_load_test);

	for (std::map<std::string, std::vector<pts_2d_conf>>::iterator it = pts_src_optimization.begin(); it != pts_src_optimization.end(); it++) {
		std::cout << "Data loaded for camera " << it->first << ": " << std::endl;
		print_vec_pts_2d_conf(it->second);
	}
	std::map<std::string, std::vector<pts_2d_conf>> pts_dst_optimization;
	multiview_optimization_multiple_points(pts_src_optimization, camera_cluster, pts_dst_optimization, consider_dist_test, num_ransac_test);
	for (std::map<std::string, std::vector<pts_2d_conf>>::iterator it = pts_dst_optimization.begin(); it != pts_dst_optimization.end(); it++) {
		std::cout << "Result for camera " << it->first << ": " << std::endl;
		print_vec_pts_2d_conf(it->second);
	}
	save_points_with_conf_multiview("image_original_conf0", 1, pts_dst_optimization);
	std::cout << "Please check if they are reasonably good." << std::endl;


	// test epipolar_2d_lines_from_anchor_point, get_3d_ray, undistort_single_point
	//std::cout << std::endl;
	//std::cout << "Testing epipolar_2d_lines_from_anchor_point, undistort_single_point............." << std::endl;
	//std::vector<pts_2d_conf> pts_src_epipolar;
	//char folder_path[1024] = "D:\\oculus\\NewAnnotationTool\\AnnotationTool_Oculus\\camera_propagation\\test\\cam330014";
	//char pts_path[1024];
	//std::sprintf(pts_path, "%s\\%s", folder_path, "00000.pose");
	//load_points_with_conf(pts_path, pts_src_epipolar, resize_factor, 1);	// load one for anchor point
	//print_vec_pts_2d_conf(pts_src_epipolar);
	//char image_path[1024];
	//std::sprintf(image_path, "%s\\%s", folder_path, "image0000.png");
	//std::cout << image_path << std::endl;
	////cv::Mat image_anchor = cv::imread(image_path, CV_LOAD_IMAGE_COLOR);   // Read the image
	////print_mat_info(image_anchor);
	////ASSERT_WITH_MSG(image_anchor.data, "Could not open or find the image");                         // Check for invalid input
	////cv::namedWindow("Anchor Image", cv::WINDOW_NORMAL);										// Create a window for display.
	////cv::circle(image_anchor, cv::Point(pts_src_epipolar[0].x, pts_src_epipolar[0].y), 10, cv::Scalar(110, 220, 0), 5);
	////cv::imshow("Anchor Image", image_anchor);												// Show our image inside it.
	////cv::resizeWindow("Anchor Image", 512, 384);
	////cv::waitKey(0);
	//mycamera camera_epipolar;
	//for (int i = 0; i < camera_cluster.size(); i++) {
	//	if (strcmp(camera_cluster[i].name.c_str(), "330014") == 0)
	//		camera_epipolar = camera_cluster[i];
	//}
	//camera_epipolar.print();
	//std::map<std::string, std::vector<double>> epipolar_lines;
	//epipolar_2d_lines_from_anchor_point(pts_src_epipolar[0], camera_epipolar, camera_cluster, epipolar_lines, consider_dist);

	//// display
	//for (std::map<std::string, std::vector<double>>::iterator it = epipolar_lines.begin(); it != epipolar_lines.end(); it++) {
	//	//std::cout << "Result for camera " << it->first << ": " << std::endl;
	//	char fileName[1024];
	//	std::sprintf(fileName, "test\\cam%s\\image0000.png", it->first.c_str());
	//	cv::Mat image_tmp = cv::imread(fileName, CV_LOAD_IMAGE_COLOR);			 // Read the image
	//	char windowName[1024];
	//	std::sprintf(windowName, "cam%s", it->first.c_str());
	//	cv::namedWindow(windowName, cv::WINDOW_NORMAL);										// Create a window for display.
	//	cv::line(image_tmp, cv::Point(0, get_y_from_2d_line(it->second, 0)), cv::Point(5120, get_y_from_2d_line(it->second, 5120)), cv::Scalar(110, 220, 0), 5);
	//	cv::imshow(windowName, image_tmp);												// Show our image inside it.
	//	cv::resizeWindow(windowName, 512, 384);
	//	cv::waitKey(0);
	//}
	//std::cout << "Please check if they are reasonably good." << std::endl;




	std::cout << std::endl;
	std::cout << "Testing done! Everything is good if you think.............." << std::endl;
	system("pause");
}

