// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu


#include "math_functions.h"
#include "type_conversion.h"
#include "debug_tool.h"

#include <opencv2/highgui.hpp> 

int main(int argc, char* argv[]) {
	// test mat
	std::vector<double> test_a = { 1,2,3 };
	//std::vector<double> test_a = { 1,2,3 };
	cv::Mat test_mat = cv::Mat(1, 2, CV_64FC3);

	std::cout << test_mat.at<cv::Vec3d>(0, 0)[0] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 0)[1] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 0)[2] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[0] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[1] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[2] << std::endl;
	test_mat.at<cv::Vec3d>(0, 0)[0] = 1;
	test_mat.at<cv::Vec3d>(0, 0)[1] = 2;
	test_mat.at<cv::Vec3d>(0, 0)[2] = 3;
	test_mat.at<cv::Vec3d>(0, 1)[0] = 4;
	test_mat.at<cv::Vec3d>(0, 1)[1] = 5;
	test_mat.at<cv::Vec3d>(0, 1)[2] = 6;
	std::cout << test_mat.at<cv::Vec3d>(0, 0)[0] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 0)[1] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 0)[2] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[0] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[1] << std::endl;
	std::cout << test_mat.at<cv::Vec3d>(0, 1)[2] << std::endl;
	//test_mat.at<double>(0, 0) = 1;
	//test_mat.at<double>(0, 0, 0) = 1;
	//test_mat.at<double>(0, 1, 0) = 2;
	//test_mat.at<double>(0, 2, 0) = 3;
	//test_mat.at<double>(0, 3, 0) = 4;
	//test_mat.at<double>(1, 0, 0) = 5;
	//test_mat.at<double>(1, 1, 0) = 6;
	//test_mat.at<double>(1, 2, 0) = 7;
	//test_mat.at<double>(1, 3, 0) = 8;
	//test_mat.at<double>(0, 0, 1) = 1;
	//test_mat.at<double>(0, 1, 1) = 2;
	//test_mat.at<double>(0, 2, 1) = 3;
	//test_mat.at<double>(0, 3, 1) = 4;
	//test_mat.at<double>(1, 0, 1) = 5;
	//test_mat.at<double>(1, 1, 1) = 6;
	//test_mat.at<double>(1, 2, 1) = 7;
	//test_mat.at<double>(1, 3, 1) = 8;
	print_mat_info(test_mat);

	// test get_median
	std::cout << std::endl;
	std::cout << "Testing get_median............." << std::endl;
	std::vector<double> test_list = { 5,3,7,2,8,1,4,0,9,6 };
	double median = get_median(test_list);
	std::cout << median << std::endl;
	ASSERT_WITH_MSG(CHECK_SCALAR_EQ(median, 4.5), "The get_median is wrong.");

	// test read_matrix
	std::cout << std::endl;
	std::cout << "Testing read_matrix............." << std::endl;
	char filename[1024];
	std::sprintf(filename, "test\\matrix.txt");	// remove surface without keypoint by cutting a plane
	FILE *califp = fopen(filename, "r");
	cv::Mat matrix(3, 3, CV_32F);
	read_matrix(califp, matrix);		// read calibration file
	print_mat(matrix);

	// test point_triangle_test_3d
	std::cout << std::endl;
	std::cout << "Testing point_triangle_test_3d............." << std::endl;
	std::vector<double> pts1 = { 1, 1, 0 };
	std::vector<double> pts2 = { 1, -1, 0 };
	std::vector<double> pts3 = { 1, 0, 1 };
	std::vector<double> pts = { 1, 0.5, 0.5 };
	if (point_triangle_test_3d(pts, pts1, pts2, pts3))
		std::cout << "The point is inside the triangle. Correct!" << std::endl;
	else
		std::cout << "The point is not inside the triangle. Wrong!" << std::endl;

	// test random_sample function
	std::cout << std::endl;
	std::cout << "Testing random_sample............." << std::endl;
	std::vector<double> samples{ 0, 0 , 1 };
	std::vector<int> selected_id_random = random_sample(samples, 100, std::time(NULL));
	print_vec(selected_id_random);

	// test generate_weighted_randomization	function
	std::cout << std::endl;
	std::cout << "Testing generate_weighted_randomization............." << std::endl;
	std::vector<double> weights{ 0, 2 , 1., 0 };
	std::vector<int> selected_id_weighted(2);
	bool all_zeros;
	for (int i = 0; i < 20; i++) {
		all_zeros = generate_weighted_randomization(weights, selected_id_weighted, std::rand());
		print_vec(selected_id_weighted);
		std::cout << "All weights are zero: " << all_zeros << std::endl;
	}

	// test l2_norm
	std::cout << std::endl;
	std::cout << "Testing l2_norm............." << std::endl;
	std::vector<double> vec_double = { 3.0, 4.0 };
	std::vector<float> vec_float = { 3.0, 4.0 };
	std::cout << l2_norm(vec_double) << std::endl;
	std::cout << l2_norm(vec_float) << std::endl;
	ASSERT_WITH_MSG(std::abs(l2_norm(vec_double) - 5) < EPS_SMALL, "The l2 norm function is wrong!");
	ASSERT_WITH_MSG(std::abs(l2_norm(vec_float) - 5) < EPS_SMALL, "The l2 norm function is wrong!");

	// test cross
	std::cout << std::endl;
	std::cout << "Testing cross............." << std::endl;
	std::vector<double> vec_double_cross1 = { 2.0, 3.0, 4.0 };
	std::vector<double> vec_double_cross2 = { 5.0, 6.0, 7.0 };
	std::vector<float> vec_float_cross1 = { 2.0, 3.0, 4.0 };
	std::vector<float> vec_float_cross2 = { 5.0, 6.0, 7.0 };
	std::vector<double> cross_double = cross(vec_double_cross1, vec_double_cross2);
	std::vector<float> cross_float = cross(vec_float_cross1, vec_float_cross2);
	std::vector<double> gt_cross = {-3, 6, -3};
	print_vec(cross_double);
	print_vec(cross_float);
	ASSERT_WITH_MSG(CHECK_VEC_EQ(cross_double, gt_cross), "The cross product is wrong!");
	ASSERT_WITH_MSG(CHECK_VEC_EQ(cross_float, double2float_vec(gt_cross)), "The cross product is wrong!");

	// test inner
	std::cout << std::endl;
	std::cout << "Testing inner............." << std::endl;
	std::vector<double> vec_double_inner1 = { 3.0, 4.0 };
	std::vector<double> vec_double_inner2 = { 4.0, -5.0 };
	std::vector<float> vec_float_inner1 = { 4.0, -5.0 };
	std::vector<float> vec_float_inner2 = { 3.0, 4.0 };
	double inner_double = inner(vec_double_inner1, vec_double_inner2);
	float inner_float = inner(vec_float_inner1, vec_float_inner2);
	std::cout << inner_double << std::endl;
	std::cout << inner_float << std::endl;
	ASSERT_WITH_MSG(std::abs(inner_double - (-8)) < EPS_SMALL, "The inner product is wrong!");
	ASSERT_WITH_MSG(std::abs(inner_float - (-8)) < EPS_SMALL, "The inner product is wrong!");

	// test get_3d_plane
	std::cout << std::endl;
	std::cout << "Testing get_plane............." << std::endl;
	pcl::PointXYZ pts_plane1(1, 0, 0);
	pcl::PointXYZ pts_plane2(0, 1, 0);
	pcl::PointXYZ pts_plane3(0, 0, 1);
	std::vector<double> plane;
	get_3d_plane(pts_plane1, pts_plane2, pts_plane3, plane);
	std::vector<double> gt_plane = { 1, 1, 1, -1 };
	print_vec(plane);
	ASSERT_WITH_MSG(CHECK_VEC_EQ(normalize_line_plane(plane), normalize_line_plane(gt_plane)), "The get_plane is wrong!");

	// test get_intersection_pts_from_lines
	std::cout << std::endl;
	std::cout << "Testing get_intersection_pts_from_lines............." << std::endl;
	std::vector<double> line_intersect1 = { 2, -1, -3 };
	std::vector<double> line_intersect2 = { 3, -1, -2 };
	cv::Point2d intersect;
	cv::Point2d gt_intersect(-1, -5);
	get_intersection_pts_from_2d_lines(line_intersect1, line_intersect2, intersect);
	print_pts2d(gt_intersect);
	ASSERT_WITH_MSG(CHECK_CV_PTS_EQ(intersect, gt_intersect), "The get_intersection_pts_from_lines is wrong!");

	// test get_projected_pts_on_line
	std::cout << std::endl;
	std::cout << "Testing get_projected_pts_on_line............." << std::endl;
	cv::Point2d pts_src_projection = { 0, 2 };
	std::vector<double> line_projection = { 1, -1, 0 };
	cv::Point2d pts_projection;
	cv::Point2d gt_pts_projection(1, 1);
	get_projected_pts_on_2d_line(pts_src_projection, line_projection, pts_projection);
	print_pts2d(pts_projection);
	ASSERT_WITH_MSG(CHECK_CV_PTS_EQ(pts_projection, gt_pts_projection), "The get_projected_pts_on_line is wrong!");


	// test get_2d_line
	std::cout << std::endl;
	std::cout << "Testing get_line............." << std::endl;
	cv::Point2d pts_line1(2, 1);
	cv::Point2d pts_line2(-1, -2);
	std::vector<double> line_vec;
	get_2d_line(pts_line1, pts_line2, line_vec);
	cv::Mat line_mat = cv::Mat(3, 1, CV_64F);
	std::vector<double> gt_line = { 1, -1, -1 };
	get_2d_line(cv2pcl_pts2d(pts_line1), cv2pcl_pts2d(pts_line2), line_mat);
	print_vec(line_vec);
	print_mat(line_mat);
	ASSERT_WITH_MSG(CHECK_VEC_EQ(normalize_line_plane(line_vec), normalize_line_plane(gt_line)), "The get_line is wrong!");
	ASSERT_WITH_MSG(CHECK_MAT_EQ(line_mat, cv::Mat(normalize_line_plane(gt_line))), "The get_line is wrong!");

	// test get_x_from_2d_line and get_y_from_2d_line
	std::cout << std::endl;
	std::cout << "Testing get_x_from_2d_line and get_y_from_2d_line............." << std::endl;
	double y1 = get_y_from_2d_line(line_vec, 2);
	double y2 = get_y_from_2d_line(line_vec, -1);
	double x1 = get_x_from_2d_line(line_vec, 1);
	double x2 = get_x_from_2d_line(line_vec, -2);
	cv::Point2d pts_got1(x1, y1);
	cv::Point2d pts_got2(x2, y2);
	print_pts2d(pts_got1);
	print_pts2d(pts_got2);
	ASSERT_WITH_MSG(CHECK_CV_PTS_EQ(pts_line1, pts_got1), "The get_x_from_2d_line or get_y_from_2d_line is wrong!");
	ASSERT_WITH_MSG(CHECK_CV_PTS_EQ(pts_line2, pts_got2), "The get_x_from_2d_line or get_y_from_2d_line is wrong!");



	std::cout << std::endl;
	std::cout << "Testing done! Everything is good............." << std::endl;
	system("pause");
}


