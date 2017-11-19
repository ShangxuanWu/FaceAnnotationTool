// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu

#include "pts_2d_conf.h"
#include "debug_tool.h"

int main(int argc, char* argv[]) {
	// ********************************************************************
	// test convert_to_point2d
	std::cout << std::endl;
	std::cout << "Testing convert_to_point2d............." << std::endl;
	pts_2d_conf pts_src(2, 3, 0.5);
	cv::Point2d pts_dst = pts_src.convert_to_point2d();
	print_pts2d(pts_dst);

	// test convert_to_pts_vec
	std::cout << std::endl;
	std::cout << "Testing convert_to_pts_vec............." << std::endl;
	std::vector<double> pts_vec1 = pts_src.convert_to_pts_vec();
	print_vec(pts_vec1);

	// test convert_to_pts_vec_conf
	std::cout << std::endl;
	std::cout << "Testing convert_to_pts_vec_conf............." << std::endl;
	std::vector<double> pts_vec2 = pts_src.convert_to_pts_vec_conf();
	print_vec(pts_vec2);

	// test print
	std::cout << std::endl;
	std::cout << "Testing print............." << std::endl;
	pts_src.print();

	std::cout << std::endl;
	std::cout << "Testing done! Everything is good............." << std::endl;
	system("pause");

	return 0;
}