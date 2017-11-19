// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu

// standard c++ library


// self-contained library
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/camera_geometry.h"
#include "../camera_propagation/IdDataParser.h"
#include "../camera_propagation/io_point.h"
#include "../camera_propagation/debug_print.h"
#include "../camera_propagation/pts_2d_conf.h"

// in-project library
#include "myheader.h"
#include "MyMesh.h"
#include "math_functions.h"
#include "pts_on_mesh.h"
#include "mesh_visualization.h"

// TODO: calculate the error in 2d and 3d

int main(int argc, char* argv[]) {
	//if (argc != 8) {
	//	fprintf(stderr, "number of arguments is %d", argc);
	//	fprintf(stderr, "calibration pose_dir mesh num_parts frame pose3dout_dir resize_factor\n");
	//	return -1;
	//}

	// define common variable
	char cmd[1024];	// store the command line
	int ret;		// store the return from running command line

					// define the input

	std::map<std::string, std::vector<pcl::PointXY>> parts2d_all_cams;
	fprintf(stderr, "reading stuff...\n");

	std::vector<mycamera> camera_cluster;
	camera_cluster.clear();
	LoadIdCalibration(argv[1], camera_cluster, false, false);	// read calibration file without skew parameter in the intrinsic
	cout << "reading calibration file successfully" << endl;

	int num_parts;
	int frame;
	double resize_factor;
	sscanf(argv[4], "%d", &num_parts);	// 51 keypoints

	sscanf(argv[5], "%d", &frame);		// frame number
	sscanf(argv[7], "%lf", &resize_factor);		// downsample scalar

	std::cout << "current path is: ";
	ret = system("cd");
	std::cout << std::endl;
	ASSERT_WITH_MSG(ret == 0, "commad line executed failed");
	std::sprintf(cmd, "if not exist %s\\del mkdir %s\\del", argv[6], argv[6]);	// remove surface without keypoint by cutting a plane
	ret = system(cmd);
	ASSERT_WITH_MSG(ret == 0, "create del folder failed");
	std::sprintf(cmd, "if not exist %s\\pose3d mkdir %s\\pose3d", argv[6], argv[6]);	// create folder for keypoints in 3D
	ret = system(cmd);
	ASSERT_WITH_MSG(ret == 0, "create pose3d folder failed");

	std::vector<std::map<std::string, pts_2d_conf>> pts_src(num_parts, std::map<std::string, pts_2d_conf>());	// [number of keypoints, map: [camera_id, 2d points (with confidence)]]
	std::map<std::string, std::vector<pts_2d_conf>> pts_load;

	char pose[1024];
	// load data from folder
	for (int i = 0; i < camera_cluster.size(); i++) {
		std::sprintf(pose, "%s\\camera%s\\%05d.pose", argv[2], camera_cluster[i].name.c_str(), frame);
		int now_part = 0;

		//std::vector<pcl::PointXY> parts2d;
		std::vector<pts_2d_conf> pts_tmp;
		load_points_with_conf(pose, pts_tmp, 10);
		ASSERT_WITH_MSG(pts_tmp.size() == num_parts, "The size of points loaded is not equal to the expected number!");
		pts_load[camera_cluster[i].name] = pts_tmp;
	}

	// convert points storage type
	std::map<std::string, pts_2d_conf> pts_multiview;
	for (int i = 0; i < pts_src.size(); i++) {
		pts_multiview.clear();
		for (std::map<std::string, std::vector<pts_2d_conf>>::iterator it = pts_load.begin(); it != pts_load.end(); it++) {
			pts_multiview[it->first] = (it->second)[i];
		}
		//std::cout << pts_multiview.size() << std::endl;
		pts_src[i] = pts_multiview;
		//std::cout << pts_src[i].size() << std::endl;
		//std::cout << pts_load.size() << std::endl;
		ASSERT_WITH_MSG(pts_load.size() == pts_src[i].size(), "The size of camera should be equal.");
	}

	// read mesh data
	std::cout << "reading mesh file" << std::endl;
	const clock_t begin_time = clock();
	MyMesh mesh(argv[3], 1);			// input reconstructed mesh
	std::cout << "It spend " << float(clock() - begin_time) / CLOCKS_PER_SEC << " second to read the mesh file" << std::endl;

	// back projection to 3d and optimization
	std::cout << "starting back projection and optimization" << std::endl << std::endl;
	std::vector<pts_on_mesh*> pts_back_3d = mesh.pts_back_projection_multiview(pts_src, camera_cluster, false);		// doesn't consider distortion
	std::vector<cv::Point3d> pts_3d;
	for (int i = 0; i < pts_back_3d.size(); i++) {
		pts_3d.push_back(pts_back_3d[i]->convert_to_point3d());
	}
	std::map<std::string, std::vector<cv::Point2d>> pts_2d_dst;
	std::cout << "startion projection from 3d to 2d" << std::endl << std::endl;
	multi_view_projection(pts_3d, camera_cluster, pts_2d_dst, false);		// projection to 2d

																			// save 2d point data
	for (int i = 0; i < camera_cluster.size(); i++) {
		std::sprintf(cmd, "if not exist testfortriangulation\\projections\\camera%s mkdir testfortriangulation\\projections\\camera%s", camera_cluster[i].name.c_str(), camera_cluster[i].name.c_str());
		ret = system(cmd);
		ASSERT_WITH_MSG(ret == 0, "camera folder in projections created failed");

		std::sprintf(pose, "testfortriangulation\\projections\\camera%s\\%010d.pose", camera_cluster[i].name.c_str(), frame);
		save_point(pts_2d_dst[camera_cluster[i].name], pose);
	}

	return 0;
}
