// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu

// standard c++ library


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

//bool inside_mouth_eye(cv::Mat& M, pcl::PointXYZ& p, std::vector<pcl::PointXY>& m) {
//	return (inside_polygon(M, p, m, 43, 51) || inside_polygon(M, p, m, 19, 25) || inside_polygon(M, p, m, 25, 31)); //mouth, right eye, left eye
//}

// implement distortion enabled

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
										
	std::vector<mycamera> camera_cluster;
	camera_cluster.clear();
	LoadIdCalibration(argv[1], camera_cluster, false);	// read calibration file
	cout << "reading calibration file successfully" << endl;
	//std::cout << camera_cluster.begin()->second.getProjectionMatrix().type() << std::endl;
	//print_mat(camera_cluster.begin()->second.getProjectionMatrix());

	int num_parts = 2;
	int frame;
	double resize_factor;
	//sscanf(argv[4], "%d", &num_parts);	// 51 keypoints

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

	////care = [28 32 34 36 37 40 43 46 49 52 55];
	////vector<int> care = {32, 34, 36, 37, 40, 43, 46, 49, 52, 55};
	//vector<int> care;
	//map<int, int> care_map;
	//for (unsigned int i = 0; i < care.size(); i++) {
	//	care_map[care[i] - 1] = 1;
	//}

	//std::vector<std::vector<std::vector<double>>> pts_estimate(num_parts, std::vector<std::vector<double>>(3, std::vector<double>()));		// save for 3d points for voting, [numKeypoints, 3, numCameras]
	std::vector<std::vector<double>> pts_visualization;													// save for 3d points for visualization, [numKeypoints, 3]
	std::vector<double> weights(num_parts, 0);
	char pose[1024];
	std::vector<std::vector<double>> rays;
	std::vector<double> ray_tmp;
	std::vector<double> camera_center;
	//int count = 0;
	//int camera_id = 330006;
	//std::vector<cv::Point3d> pts_src;				
	//for (; count < 1; count++) {
		//cv::Mat i = Ms.find(camera_id)->second;
		//int id = Ms.find(camera_id)->first;
	//for (map<int, cv::Mat>::iterator i = Ms.begin(); count < 1; i++) {
	for (std::map<std::string, mycamera>::iterator i = camera_cluster.begin(); i != camera_cluster.end(); i++) {
		print_mat(i->second.getProjectionMatrix());

		std::sprintf(pose, "%s\\camera%s\\%05d.pose", argv[2], i->first.c_str(), frame);
		fprintf(stderr, "processing frame 1, the 2d keypoint file is in the %s\n", pose);
		cout << "reading 2d keypoint file" << endl;
		FILE* in = fopen(pose, "r");
		if (in == NULL) {
			fprintf(stderr, "skipping %s\n", i->first.c_str());
			continue;
		}
		double x, y, conf;
		int now_part = 0;

		std::vector<pcl::PointXY> parts2d;
		while (fscanf(in, "%lf %lf %lf\n", &x, &y, &conf) > 0) {	// read all keypoint coordinates and confidence
			fprintf(stderr, "\npart %d\n", now_part);
			x *= resize_factor;		// change to original scale
			y *= resize_factor;
			pcl::PointXY p;
			p.x = x;
			p.y = y;
			cout << "2D coordinate in the original image is (" << x << ", " << y << ")" << endl;
			parts2d.push_back(p);
			/*if( care_map.find(now_part) == care_map.end() || conf < 0.9 ){
			now_part++;
			continue;
			}*/

			// TODO: check the correctness of this
			ray_tmp.clear();
			camera_center.clear(); 
			pts_on_mesh* pom = mesh.project_prediction(x, y, i->second.getProjectionMatrix(), conf, ray_tmp, camera_center);		// retrieve the 3d point coordinate
			ASSERT_WITH_MSG(ray_tmp.size() == 4, "The size of the output ray is not correct. Please check!");
			ASSERT_WITH_MSG(camera_center.size() == 4, "The size of the camera center is not correct. Please check!");
			rays.push_back(ray_tmp);

			//cv::Mat center_mat(4, 1, CV_64F);
			//center_mat = cv::Mat(camera_center);
			//std::cout << "camera center is " << std::endl;
			//print_vec(camera_center);
			////print_mat(center_mat);
			//std::cout << std::endl;

			//std::cout << "projection matrix is " << std::endl;
			//cv::Mat projection_matrix(3, 4, CV_64F);
			//projection_matrix = i->second.getProjectionMatrix();
			//print_mat(projection_matrix);
			//std::cout << std::endl;

			//cv::Mat result(3, 1, CV_64F);
			//result = projection_matrix * center_mat;
			//std::cout << "camera center in 2d is " << std::endl;
			//print_mat(result);
			//std::cout << std::endl;
			//std::vector<double> pts_2d;
			//for (int k = 0; k < result.rows; k++) {
			//	pts_2d.push_back(result.at<double>(k, 0));
			//}
			//std::cout << "camera center in 2d is " << std::endl;
			//print_vec(pts_2d);
			//std::cout << std::endl;

			//std::cout << l2_norm(pts_2d) << std::endl;

			//ASSERT_WITH_MSG(abs(l2_norm(pts_2d)) < EPS, "The output of camera center is not correct!");


			/*if( pom.conf < 0.9 ){
			now_part++;
			continue;
			}*/
			//fprintf(stderr, "%d %g %g %g %g\n", acc, pom.x, pom.y, pom.z, conf);
			cout << "3D coordinate in the original image is (" << pom->x << ", " << pom->y << ", " << pom->z << ")" << endl;
			//pts_estimate[now_part][0].push_back(pom->x);
			//pts_estimate[now_part][1].push_back(pom->y);
			//pts_estimate[now_part][2].push_back(pom->z);
			//weights[now_part] += pom->conf;
			//pts_src.push_back(cv::Point3d(pom.x, pom.y, pom.z));	// save 3d point to the vector
			
			// for visualization
			std::vector<double> single_point;
			single_point.push_back(pom->x);
			single_point.push_back(pom->y);
			single_point.push_back(pom->z);
			pts_visualization.push_back(single_point);

			//print_vec(single_point);
			//print_vec(pts_visualization[0]);
			now_part++;
			if (now_part == 2)
				break;
		}
		std::fclose(in);
		std::sprintf(cmd, "if not exist %s\\%s mkdir %s\\%s", argv[6], i->first.c_str(), argv[6], i->first.c_str());
		ret = system(cmd);
		ASSERT_WITH_MSG(ret == 0, "camera folder in wrapped created failed");
		parts2d_all_cams[i->first] = parts2d;
	}
	std::vector<double> origin;
	origin.push_back(camera_center[0]);
	origin.push_back(camera_center[1]);
	origin.push_back(camera_center[2]);
	pts_visualization.push_back(origin);
	//print_vec(pts_visualization[0]);

	//test git
	// test visualization
	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr line_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	//boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	//get_cloud_from_points(pts_visualization, keypoints_cloud_ptr);
	//pcl::PointXYZ pts_start_tmp(camera_center[0], camera_center[1], camera_center[2]);
	//std::vector<pcl::PointXYZ> pts_start;
	//std::vector<uint32_t> range;
	//for (int i = 0; i < num_parts; i++) {
	//	pts_start.push_back(pts_start_tmp);
	//	range.push_back(10000);
	//}
	//get_cloud_from_lines(rays, line_cloud_ptr, pts_start, range);
	//viewer = keypoint_line_mesh_visualization(keypoints_cloud_ptr, line_cloud_ptr, mesh.cloud);
	//while (!viewer->wasStopped())
	//{
	//	viewer->spinOnce(100);
	//	boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	//}

	// save projected 2d points to all cameras
	std::map<std::string, std::vector<cv::Point2d>> pts_dst;
	project(pts_src, pts_dst, camera_cluster);	
	for (std::map<std::string, mycamera>::iterator camera_temp = camera_cluster.begin(); camera_temp != camera_cluster.end(); camera_temp++) {
		std::sprintf(cmd, "if not exist testfortriangulation\\projections\\camera%s mkdir testfortriangulation\\projections\\camera%s", camera_temp->second.name.c_str(), camera_temp->second.name.c_str());
		ret = system(cmd);
		ASSERT_WITH_MSG(ret == 0, "camera folder in projections created failed");

		std::sprintf(pose, "testfortriangulation\\projections\\camera%s\\%05d.pose", camera_temp->second.name.c_str(), frame);
		save_point(pts_dst[camera_temp->second.name], pose);
	}


	////find closest triangle and project...
	//std::vector<int> closest;
	//sprintf(cmd, "%s\\pose3d\\%05d.pose3d", argv[6], frame);
	//FILE *out = fopen(cmd, "w");
	//ASSERT_WITH_MSG(out != NULL, "pose3d file cannot be opened for saving 3d coordinate");
	//int neighbors = 1;
	//std::vector<int> pointIdxNKNSearch(neighbors);
	//std::vector<float> pointNKNSquaredDistance(neighbors);
	//pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	//kdtree.setInputCloud(mesh.cloud);
	//std::vector<pcl::PointXYZ> svd_m;
	//for (int i = 0; i < num_parts; i++) {

	//	if (pts_estimate[i][0].size() == 0) {
	//		fprintf(out, "-1 -1 -1 -1 -1 -1\n");
	//		closest.push_back(-1);
	//		continue;
	//	}
	//	/*
	//	for(int j = 0; j < 3; j++){
	//	pts_estimate[i][j] /= weights[i];
	//	}
	//	*/
	//	pcl::PointXYZ pts;
	//	mean_shift(pts_estimate[i], pts);
	//	//pts.x = get_median(pts_estimate[i][0]);
	//	//pts.y = get_median(pts_estimate[i][1]);
	//	//pts.z = get_median(pts_estimate[i][2]);
	//	kdtree.nearestKSearch(pts, neighbors, pointIdxNKNSearch, pointNKNSquaredDistance);
	//	int closest_id = pointIdxNKNSearch[0];
	//	closest.push_back(closest_id);
	//	//fprintf(stderr, "k %g %g %g %g\n", pts.x, pts.y, pts.z, sqrt(pointNKNSquaredDistance[0]));
	//	//fprintf(stderr, "p %g %g %g %g\n", mesh.cloud->points[closest_id].x, mesh.cloud->points[closest_id].y, mesh.cloud->points[closest_id].z, weights[i]);
	//	//fprintf(out, "%g %g %g %g %d %d\n", mesh.cloud->points[closest_id].x, mesh.cloud->points[closest_id].y, mesh.cloud->points[closest_id].z, weights[i], closest_id, mesh.first_plane_of_vertices[closest_id]);
	//	fprintf(out, "%g %g %g %g\n", pts.x, pts.y, pts.z, weights[i]);
	//	svd_m.push_back(pts);
	//}
	//fclose(out);

	//Eigen::MatrixXf mymatrix(svd_m.size(), 4);
	//for (int i = 0; i < svd_m.size(); i++) {
	//	mymatrix(i, 0) = svd_m[i].x;
	//	mymatrix(i, 1) = svd_m[i].y;
	//	mymatrix(i, 2) = svd_m[i].z;
	//	mymatrix(i, 3) = 1;
	//}
	//Eigen::JacobiSVD<Eigen::MatrixXf> svd(mymatrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
	//std::vector<float> cutting_plane;
	//auto matV = svd.matrixV();
	//for (int i = 0; i < 4; i++) {
	//	cutting_plane.push_back(matV(i, 3)); //last column
	//}

	//float acc = 0;
	//acc += svd_m[0].x * cutting_plane[0]; //svd_m[0] is left eyebrow
	//acc += svd_m[0].y * cutting_plane[1]; //svd_m[0] is left eyebrow
	//acc += svd_m[0].z * cutting_plane[2]; //svd_m[0] is left eyebrow
	//acc += cutting_plane[3];
	//cutting_plane[3] -= acc; //if the plane was distance d away, now make it 0 away
	//acc = 0;
	//acc += svd_m[13].x * cutting_plane[0]; //svd_m[13] is tip of nose
	//acc += svd_m[13].y * cutting_plane[1];
	//acc += svd_m[13].z * cutting_plane[2];
	//acc += cutting_plane[3];



	//sprintf(cmd, "%s\\del\\%05d.obj", argv[6], frame);
	//out = fopen(cmd, "w");
	//ASSERT_WITH_MSG(out != NULL, "obj file cannot be opened for saving warppred mesh");
	//unordered_set<int> faces_to_delete;
	//int pts = mesh.cloud->size();
	//map<int, int> remap_vertices;
	//for (int i = 0; i < pts; i++) {
	//	float acc2 = 0;
	//	acc2 += cutting_plane[0] * mesh.cloud->points[i].x;
	//	acc2 += cutting_plane[1] * mesh.cloud->points[i].y;
	//	acc2 += cutting_plane[2] * mesh.cloud->points[i].z;
	//	acc2 += cutting_plane[3];
	//	if (acc * acc2 < 0 || inside_mouth_eye(Ms[330030], mesh.cloud->points[i], parts2d_all_cams["330030"])) { // not on same side as the detected points, throw away
	//		for (int j = 0; j < mesh.planes_of_vertices[i].size(); j++) {
	//			faces_to_delete.insert(mesh.planes_of_vertices[i][j]);
	//		}
	//	}
	//	else {
	//		remap_vertices[i] = remap_vertices.size();
	//		fprintf(out, "v %g %g %g\n", mesh.cloud->points[i].x, mesh.cloud->points[i].y, mesh.cloud->points[i].z);
	//	}
	//}

	//for (int i = 0; i < mesh.plane_pts_idx.size(); i++) {
	//	if (faces_to_delete.find(i) == faces_to_delete.end()) {
	//		int ok = 1;
	//		for (int j = 0; j < 3; j++) {
	//			if (remap_vertices.find(mesh.plane_pts_idx[i][j]) == remap_vertices.end()) {
	//				ok = 0;
	//				break;
	//			}
	//		}
	//		if (ok == 1)
	//			fprintf(out, "f %d %d %d\n", remap_vertices[mesh.plane_pts_idx[i][0]], remap_vertices[mesh.plane_pts_idx[i][1]], remap_vertices[mesh.plane_pts_idx[i][2]]);
	//	}
	//}
	//fclose(out);

	////for (map<int, Mat>::iterator j = Ms.begin(); j != Ms.end(); j++) {
	//count = 0;
	//for (map<int, Mat>::iterator j = Ms.begin(); count < 1; j++, count++) {
	//	sprintf(cmd, "%s\\%d\\%05d.pose", argv[6], j->first, frame);
	//	FILE *out = fopen(cmd, "w");
	//	assert(out != NULL);
	//	for (int i = 0; i < num_parts; i++) {
	//		if (closest[i] == -1) {
	//			fprintf(out, "-1 -1 -1\n");
	//			continue;
	//		}
	//		Mat P(4, 1, CV_32F);
	//		P.at<float>(0, 0) = mesh.cloud->points[closest[i]].x;
	//		P.at<float>(1, 0) = mesh.cloud->points[closest[i]].y;
	//		P.at<float>(2, 0) = mesh.cloud->points[closest[i]].z;
	//		P.at<float>(3, 0) = 1;
	//		P = j->second * P;
	//		fprintf(out, "%g %g %g\n", P.at<float>(0, 0) / P.at<float>(2, 0) / resize_factor, P.at<float>(1, 0) / P.at<float>(2, 0) / resize_factor, weights[i]);
	//	}
	//	fclose(out);
	//}
	return 0;
}
