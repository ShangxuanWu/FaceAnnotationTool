// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu

#include "io_point.h"
#include "mycamera.h"
#include "camera_geometry.h"
#include "IdDataParser.h"
#include "pts_2d_conf.h"

#define consider_dist	true
#define consider_skew	false
#define resize_factor	1
#define number_pts_load	3

int main(int argc, char* argv[]) {
	// ********************************************************************
	// read calibration file
	std::vector<mycamera> camera_cluster;
	camera_cluster.clear();
	LoadIdCalibration("calibration.txt", camera_cluster, consider_skew, consider_dist);		// consider distortion	
	
	// test save_points_with_conf and save_points_with_conf_multiview
	int frame = 0;
	std::map<std::string, std::vector<pts_2d_conf>> pts_src_load;
	load_points_with_conf_multiview("test", camera_cluster, frame, pts_src_load, resize_factor, number_pts_load);
	save_points_with_conf_multiview("test_save", frame, pts_src_load);


	std::cout << std::endl;
	std::cout << "Testing done! Everything is good............." << std::endl;
	system("pause");

	return 0;
}