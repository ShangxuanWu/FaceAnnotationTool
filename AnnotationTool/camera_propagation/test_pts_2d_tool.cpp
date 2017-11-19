// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu

#include "pts_2d_tool.h"
#include "pts_2d_conf.h"

int main(int argc, char* argv[]) {
	// ********************************************************************
	// test convert_to_pts_2d_conf
	std::cout << std::endl;
	std::cout << "Testing convert_to_pts_2d_conf............." << std::endl;
	pts_2d_tool pts_src(2, 3, 0.5, true);
	pts_2d_conf pts_dst = pts_src.convert_to_pts_2d_conf();
	pts_dst.print();

	// test print
	std::cout << std::endl;
	std::cout << "Testing print............." << std::endl;
	pts_src.print();

	// test constructor
	std::cout << std::endl;
	std::cout << "Testing constructor............." << std::endl;
	pts_2d_tool pts_cons1(pts_dst, true);
	pts_cons1.print();
	pts_2d_tool pts_cons2(pts_dst);
	pts_cons2.print();

	std::cout << std::endl;
	std::cout << "Testing done! Everything is good............." << std::endl;
	system("pause");

	return 0;
}