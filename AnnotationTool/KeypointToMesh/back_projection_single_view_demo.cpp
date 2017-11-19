// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu

// self-contained library
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/camera_geometry.h"
#include "../camera_propagation/IdDataParser.h"
#include "../camera_propagation/io_point.h"
#include "../camera_propagation/debug_print.h"

// in-project library
#include "myheader.h"
#include "MyMesh.h"
#include "math_functions.h"
#include "pts_on_mesh.h"
#include "mesh_visualization.h"

int main(int argc, char* argv[]) {
	if (argc != 8) {
		fprintf(stderr, "number of arguments is %d", argc);
		fprintf(stderr, "calibration pose_dir mesh num_parts frame pose3dout_dir resize_factor\n");
		return -1;
	}

	char cmd[1024];	// store the command line
	int ret;		// store the return from running command line

					
	// reading parameters and data from disk
	//std::map<int, cv::Mat> Ms;
	std::map<std::string, std::vector<pcl::PointXY>> parts2d_all_cams;
	fprintf(stderr, "reading stuff...\n");
	//read_calibration(argv[1], Ms);		// read calibration file

	std::map<std::string, mycamera> camera_cluster;
	camera_cluster.clear();
	LoadIdCalibration(argv[1], camera_cluster, false);	// read calibration file
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

	std::cout << "reading mesh file" << std::endl;
	const clock_t begin_time = clock();
	MyMesh mesh(argv[3], 1);			// input reconstructed mesh
	std::cout << "It spend " << float(clock() - begin_time) / CLOCKS_PER_SEC << " second to read the mesh file" << std::endl;

	std::vector<std::vector<std::vector<double>>> pts_estimate(num_parts, std::vector<std::vector<double>>(3, std::vector<double>()));		// save for 3d points for voting, [numKeypoints, 3, numCameras]
	std::vector<std::vector<double>> pts_visualization;													// save for 3d points for visualization, [numKeypoints, 3]
	std::vector<double> weights(num_parts, 0);
	char pose[1024];
	std::vector<std::vector<double>> rays;
	std::vector<double> ray_tmp;
	std::vector<double> camera_center;

	for (std::map<std::string, mycamera>::iterator i = camera_cluster.begin(); i != camera_cluster.end(); i++) {
		print_mat(i->second.getProjectionMatrix());

		std::sprintf(pose, "%s\\camera%s\\%05d.pose", argv[2], i->first.c_str(), frame);
		fprintf(stderr, "processing frame 1, the 2d keypoint file is in the %s\n", pose);
		cout << "reading 2d keypoint file" << endl;
		FILE* in = fopen(pose, "r");
		if (in == NULL) {
			fprintf(stderr, "skipping %d\n", i->first.c_str());
			continue;
		}
		double x, y, conf;
		int now_part = 0;

		std::vector<pcl::PointXY> parts2d;
		while (fscanf(in, "%lf %lf %lf\n", &x, &y, &conf) > 0) {	// read all keypoint coordinates and confidence
			fprintf(stderr, "part %d\n", now_part);
			x *= resize_factor;		// change to original scale
			y *= resize_factor;
			pcl::PointXY p;
			p.x = x;
			p.y = y;
			cout << "2D coordinate in the original image is (" << x << ", " << y << ")" << endl;
			parts2d.push_back(p);

			ray_tmp.clear();
			camera_center.clear();
			pts_on_mesh* pom = mesh.project_prediction(x, y, i->second.getProjectionMatrix(), conf, ray_tmp, camera_center);		// retrieve the 3d point coordinate
			ASSERT_WITH_MSG(ray_tmp.size() == 4, "The size of the output ray is not correct. Please check!");
			ASSERT_WITH_MSG(camera_center.size() == 4, "The size of the camera center is not correct. Please check!");
			rays.push_back(ray_tmp);

			cout << "3D coordinate in the original image is (" << pom->x << ", " << pom->y << ", " << pom->z << ")" << endl;
			pts_estimate[now_part][0].push_back(pom->x);
			pts_estimate[now_part][1].push_back(pom->y);
			pts_estimate[now_part][2].push_back(pom->z);
			weights[now_part] += pom->conf;

			// for visualization
			std::vector<double> single_point;
			single_point.push_back(pom->x);
			single_point.push_back(pom->y);
			single_point.push_back(pom->z);
			pts_visualization.push_back(single_point);

			now_part++;
		}
		std::fclose(in);
		std::sprintf(cmd, "if not exist %s\\%d mkdir %s\\%d", argv[6], i->first.c_str(), argv[6], i->first.c_str());
		ret = system(cmd);
		ASSERT_WITH_MSG(ret == 0, "camera folder in wrapped created failed");
		parts2d_all_cams[i->first] = parts2d;
	}
	std::vector<double> origin;
	origin.push_back(camera_center[0]);
	origin.push_back(camera_center[1]);
	origin.push_back(camera_center[2]);
	pts_visualization.push_back(origin);

	// test visualization
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr line_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	get_cloud_from_points(pts_visualization, keypoints_cloud_ptr);
	pcl::PointXYZ pts_start_tmp(camera_center[0], camera_center[1], camera_center[2]);
	std::vector<pcl::PointXYZ> pts_start;
	std::vector<uint32_t> range;
	for (int i = 0; i < num_parts; i++) {
		pts_start.push_back(pts_start_tmp);
		range.push_back(10000);
	}
	get_cloud_from_lines(rays, line_cloud_ptr, pts_start, range);
	viewer = keypoint_line_mesh_visualization(keypoints_cloud_ptr, line_cloud_ptr, mesh.cloud);
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}

	return 0;
}
