// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu


// in-project library
#include "MyMesh.h"
#include "pts_on_mesh.h"


// pcl library
#include <pcl/conversions.h>
#include <pcl/PolygonMesh.h>
#include <pcl/io/auto_io.h>


// self_contained library
#include "../camera_propagation/debug_tool.h"
#include "../camera_propagation/mycamera.h"
#include "../camera_propagation/camera_geometry.h"
#include "../camera_propagation/pts_2d_conf.h"
#include "../camera_propagation/pts_3d_conf.h"
#include "../camera_propagation/type_conversion.h"
#include "../camera_propagation/math_functions.h"


MyMesh::MyMesh(char* filename, int scale) {
	pcl::PolygonMesh::Ptr mesh(new pcl::PolygonMesh);
	pcl::io::load(filename, *mesh.get());
	cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::fromPCLPointCloud2(mesh->cloud, *cloud.get());
	for (unsigned int i = 0; i < cloud->size(); i++) {
		cloud->points[i].x *= scale;		//rescale the point to the original scale
		cloud->points[i].y *= scale;
		cloud->points[i].z *= scale;
	}

	int num_vertices = cloud->size();
	int num_planes = mesh->polygons.size();
	planes_of_vertices.resize(num_vertices);
	first_plane_of_vertices.resize(num_vertices);
	planes.resize(num_planes);

	for (int i = 0; i < num_planes; i++) {
		ASSERT_WITH_MSG(mesh->polygons[i].vertices.size() == 3, "The vertices of mesh should have 3 coordinates");			// ensure mesh is in 3d
		std::vector<int> idxs(3);
		for (int j = 0; j < 3; j++) {
			idxs[j] = mesh->polygons[i].vertices[j];			// find the id of vertices corresponding to all polygon (triangle)
		}
		//std::cout << "idxs is:" << std::endl;
		//print_vec(idxs);

		plane_pts_idx.push_back(idxs);							// idxs is three points id for current plane, plane_pts_idx store the vertices id for all planes

		Eigen::Matrix3f A, B;									// A stores 3 3d points for each polygon
		A << cloud->points[idxs[0]].x, cloud->points[idxs[0]].y, cloud->points[idxs[0]].z,
			cloud->points[idxs[1]].x, cloud->points[idxs[1]].y, cloud->points[idxs[1]].z,
			cloud->points[idxs[2]].x, cloud->points[idxs[2]].y, cloud->points[idxs[2]].z;		// cloud stores the exact 3d coordiante indexed by vertices id
		B << 0, 0, 0, 1, 0, 0, 0, 1, 0;
		Eigen::ColPivHouseholderQR<Eigen::Matrix3f> dec(A);
		plane_projection.push_back(dec.solve(B));			// keep
		if (dec.rank() != 3) {
			plane_projection_good.push_back(false);
		}
		else {
			plane_projection_good.push_back(true);
		}

		double normal_length = get_3d_plane(cloud->points[idxs[0]], cloud->points[idxs[1]], cloud->points[idxs[2]], planes[i]);		// planes is the 4d vector representing that plane
		if (normal_length > EPS_SMALL) {			// if the length is very small, then the three points forming that plane is in the same straight line
			first_plane_of_vertices[idxs[0]] = i;
			for (int j = 0; j < 3; j++) {
				planes_of_vertices[idxs[j]].push_back(i);		// planes of vertices stores which plane this points belongs to indexed by points id
			}
		}
	}
}


MyMesh::~MyMesh() {
}


// TODO: test for correctness
pts_on_mesh* MyMesh::pts_back_projection_single_view(pts_2d_conf& pts_2d, mycamera& mycamera, const bool consider_dist) {
	std::vector<double> ray;
	cv::Point3d C;
	get_3d_ray(pts_2d, mycamera, C, ray, consider_dist);
	return this->get_pts_on_mesh(C, ray, pts_2d.conf);
}

// TODO: test for correctness
pts_on_mesh* MyMesh::pts_back_projection_multiview(std::map<std::string, pts_2d_conf>& pts_2d, std::vector<mycamera>& camera_cluster, const bool consider_dist) {
	ASSERT_WITH_MSG(consider_dist == false, "Back projection doesn't support undistortion right now. Please undistort first!");
	ASSERT_WITH_MSG(pts_2d.size() == camera_cluster.size(), "The size of the input points and camera is not equal.");
	double error_minimum = -1;
	int minimum_camera_id = -1;
	double error_tmp;
	for (int i = 0; i < camera_cluster.size(); i++) {
		pts_2d_conf pts_tmp = pts_2d[camera_cluster[i].name];
		//pts_tmp.print();
		pts_on_mesh* pts_3d = pts_back_projection_single_view(pts_tmp, camera_cluster[i], consider_dist);	// obtain 3d point on the mesh from single view
		//pts_3d->print();

		// calculate error
		cv::Point3d pts_src(pts_3d->x, pts_3d->y, pts_3d->z);
		//error_tmp = calculate_projection_error(pts_src, camera_cluster, pts_2d, true, consider_dist);
		error_tmp = calculate_projection_error(pts_src, camera_cluster, pts_2d, true);

		if (error_minimum == -1 || error_tmp < error_minimum) {
			error_minimum = error_tmp;
			minimum_camera_id = i;
		}
	}

	return pts_back_projection_single_view(pts_2d[camera_cluster[minimum_camera_id].name], camera_cluster[minimum_camera_id], consider_dist);
}

// TODO: test for correctness
std::vector<pts_on_mesh*> MyMesh::pts_back_projection_multiview(std::vector<std::map<std::string, pts_2d_conf>>& pts_2d, std::vector<mycamera>& camera_cluster, const bool consider_dist) {
	std::cout << "starting back projection for multiple points from multiview camera" << std::endl << std::endl;
	ASSERT_WITH_MSG(consider_dist == false, "Back projection doesn't support undistortion right now. Please undistort first!");
	std::vector<pts_on_mesh*> pts_dst;
	//std::cout << pts_2d.size() << std::endl;
	for (int i = 0; i < pts_2d.size(); i++) {
		std::cout << "working on view " << i << ", ";
		ASSERT_WITH_MSG(pts_2d[i].size() == camera_cluster.size(), "The size of the input points and camera is not equal.");
		pts_dst.push_back(pts_back_projection_multiview(pts_2d[i], camera_cluster, consider_dist));
	}
	return pts_dst;
}

// TODO: test for correctness
pts_on_mesh* MyMesh::get_pts_on_mesh(cv::Point3d C_src, std::vector<double>& ray, double conf) {
	//std::cout << "finding the closest point on the mesh given a 3d point." << std::endl;
	double closest = -1;
	std::vector<double> cpts;
	int tri_id = -1;
	int planes = this->planes.size();
	double t;
	int k;
	std::vector<double> pts_tmp;
	std::vector<double> C = cv2vec_pts3d(C_src);

	// go though all planes
	for (k = 0; k < planes; k++) {
		// plane is plane[0] * x + plane[1] * y + plane[2] * z + plane[3] = 0
		// expected 3d point is [C[0] + t * ray[0], C[1] + t * ray[1], C[2] + t * ray[2]], t is the parameter we need to find
		std::vector<double>& plane = this->planes[k];
		double n = 0, d = 0;
		for (int l = 0; l < 3; l++) {
			n += C[l] * plane[l];
			d += ray[l] * plane[l];
		}
		n += plane[3];
		if (abs(d) < EPS_MYSELF) {
			continue;
		}
		t = -n / d;
		if (t < 0) {
			continue;
		}

		pts_tmp.clear();
		for (int i = 0; i < 3; i++) {
			pts_tmp.push_back(t * ray[i] + C[i]);	// the 3d point
		}

		std::vector<double> pts_append = pts_tmp;
		pts_append.push_back(1.0);

		ASSERT_WITH_MSG(abs(inner(pts_append, plane)) < EPS_MYSELF, "Point found is not even on the extended triangle plane. The offset is " + std::to_string(abs(inner(pts_append, plane))));		// check the point with parameter t is inside the current plane

		std::vector<int> vec_id = this->plane_pts_idx[k];
		std::vector<double> tri_a, tri_b, tri_c;
		tri_a.push_back(this->cloud->points[vec_id[0]].x);
		tri_a.push_back(this->cloud->points[vec_id[0]].y);
		tri_a.push_back(this->cloud->points[vec_id[0]].z);
		tri_b.push_back(this->cloud->points[vec_id[1]].x);
		tri_b.push_back(this->cloud->points[vec_id[1]].y);
		tri_b.push_back(this->cloud->points[vec_id[1]].z);
		tri_c.push_back(this->cloud->points[vec_id[2]].x);
		tri_c.push_back(this->cloud->points[vec_id[2]].y);
		tri_c.push_back(this->cloud->points[vec_id[2]].z);

		// check if the point is inside the triangle
		if (point_triangle_test_3d(pts_tmp, tri_a, tri_b, tri_c)) {
			if (tri_id == -1 || t < closest) {		// find the exact right plane which intersect with the ray and the intersection point is inside the plane.
				closest = t;
				tri_id = k;
				cpts = pts_tmp;
			}
		}
	}

	//std::cout << "plane id is:" << tri_id << std::endl;
	//std::cout << std::endl;
	//std::cout << "plane is:" << std::endl;
	//print_vec(this->planes[tri_id]);
	//std::cout << std::endl;

	//std::cout << "id of 3 points around the plane are" << std::endl;
	std::vector<int> ids = this->plane_pts_idx[tri_id];
	//std::cout << ids[0] << ", " << ids[1] << ", " << ids[2] << std::endl;
	//std::cout << std::endl;

	//std::cout << "The 3d coordinate of three points in the intersecting plane are:" << std::endl;
	//std::vector<double> pts1, pts2, pts3;
	//pts1.push_back(this->cloud->points[ids[0]].x);
	//pts1.push_back(this->cloud->points[ids[0]].y);
	//pts1.push_back(this->cloud->points[ids[0]].z);
	//pts2.push_back(this->cloud->points[ids[1]].x);
	//pts2.push_back(this->cloud->points[ids[1]].y);
	//pts2.push_back(this->cloud->points[ids[1]].z);
	//pts3.push_back(this->cloud->points[ids[2]].x);
	//pts3.push_back(this->cloud->points[ids[2]].y);
	//pts3.push_back(this->cloud->points[ids[2]].z);
	//print_vec(pts1);
	//print_vec(pts2);
	//print_vec(pts3);
	//std::cout << std::endl;

	//std::cout << "t is:" << t << std::endl;
	//std::cout << std::endl;

	//print_vec(cpts);
	//std::cout << std::endl;

	pts_on_mesh* ptr_mesh;
	if (tri_id == -1) {
		conf = 0;
		fprintf(stderr, "no projected triangles!\n");

		ptr_mesh = new pts_on_mesh(-1, 0.0, 0.0, 0.0, conf);
		return ptr_mesh;
	}
	//ptr_mesh = new pts_on_mesh(this->plane_pts_idx[tri_id][0], cpts[0], cpts[1], cpts[2], conf);
	pts_3d_conf pts_3d(cpts[0], cpts[1], cpts[2], conf);
	//std::cout << "final 3d points is" << std::endl;
	//pts_3d.print();
	pts_on_mesh* pts_mesh = find_closest_pts_on_mesh(pts_3d, tri_id);
	pts_mesh->print();
	return pts_mesh;
}

// TODO: test for correctness
pts_on_mesh* MyMesh::find_closest_pts_on_mesh(pts_3d_conf& pts_3d) {
	//std::cout << "matching 3d point found with points on the mesh." << std::endl;
	ASSERT_WITH_MSG(pts_3d.conf >= 0 && pts_3d.conf <= 1, "The confidence should be in the range of [0, 1].");
	int closest_id = -1;
	double closest_distance = -1;
	double dist_tmp;
	for (int i = 0; i < this->cloud->points.size(); i++) {
		dist_tmp = sqrt((this->cloud->points[i].x - pts_3d.x) * (this->cloud->points[i].x - pts_3d.x) + (this->cloud->points[i].y - pts_3d.y) * (this->cloud->points[i].y - pts_3d.y) + (this->cloud->points[i].z - pts_3d.z) * (this->cloud->points[i].z - pts_3d.z));
		if (closest_distance == -1 || dist_tmp < closest_distance) {
			closest_distance = dist_tmp;
			closest_id = i;
		}
	}
	return new pts_on_mesh(closest_id, this->cloud->points[closest_id].x, this->cloud->points[closest_id].y, this->cloud->points[closest_id].z, pts_3d.conf);
}

// TODO: test for correctness
pts_on_mesh* MyMesh::find_closest_pts_on_mesh(pts_3d_conf& pts_3d, int plane_id) {
	//std::cout << "matching 3d point found with points on the mesh." << std::endl;
	ASSERT_WITH_MSG(pts_3d.conf >= 0 && pts_3d.conf <= 1, "The confidence should be in the range of [0, 1].");
	int number_vertices = this->plane_pts_idx[plane_id].size();
	//std::cout << "number vertices: " << number_vertices << std::endl;
	int vertice_id_tmp;
	int closest_id = -1;
	double closest_distance = -1;
	double dist_tmp;
	for (int i = 0; i < number_vertices; i++) {
		vertice_id_tmp = this->plane_pts_idx[plane_id][i];
		//std::cout << vertice_id_tmp << std::endl;
		dist_tmp = sqrt((this->cloud->points[vertice_id_tmp].x - pts_3d.x) * (this->cloud->points[vertice_id_tmp].x - pts_3d.x) + (this->cloud->points[vertice_id_tmp].y - pts_3d.y) * (this->cloud->points[vertice_id_tmp].y - pts_3d.y) + (this->cloud->points[vertice_id_tmp].z - pts_3d.z) * (this->cloud->points[vertice_id_tmp].z - pts_3d.z));
		//std::cout << dist_tmp << std::endl;
		if (closest_distance == -1 || dist_tmp < closest_distance) {
			closest_distance = dist_tmp;
			closest_id = this->plane_pts_idx[plane_id][i];
		}
	}
	//std::cout << closest_id << std::endl;
	return new pts_on_mesh(closest_id, this->cloud->points[closest_id].x, this->cloud->points[closest_id].y, this->cloud->points[closest_id].z, pts_3d.conf);
}



// deprecated
// optimize the multiview triangulation in 3d space, current strategy is to select one best 3d point from all view
// support only one point
// optimization involved
//pts_on_mesh* MyMesh::pts_back_projection_multiview(std::map<std::string, pts_2d_conf>& pts_2d, std::map<std::string, mycamera>& camera_cluster) {
//	ASSERT_WITH_MSG(pts_2d.size() == camera_cluster.size(), "The size of the input points and camera is not equal.");
//	double error_minimum = -1;
//	std::string minimum_camera_id;
//	double error_tmp;
//	for (std::map<std::string, mycamera>::iterator it = camera_cluster.begin(); it != camera_cluster.end(); it++) {
//		pts_2d_conf pts_tmp = pts_2d[it->first];
//		pts_on_mesh* pts_3d = pts_back_projection_single_view(pts_tmp, it->second);	// obtain 3d point on the mesh from single view
//																								
//		cv::Point3d pts_src(pts_3d->x, pts_3d->y, pts_3d->z);
//		error_tmp = calculate_projection_error(pts_src, camera_cluster, pts_2d, true);// calculate error
//
//		if (error_minimum == -1 || error_tmp < error_minimum) {
//			error_minimum = error_tmp;
//			minimum_camera_id = it->first;
//		}
//	}
//
//	return pts_back_projection_single_view(pts_2d[minimum_camera_id], camera_cluster[minimum_camera_id]);
//}

// deprecated
// multiple points from all views
// each pts_2d could be only x, y or includes confidence
// optimization involved
//std::vector<pts_on_mesh*> MyMesh::pts_back_projection_multiview(std::vector<std::map<std::string, pts_2d_conf>>& pts_2d, std::map<std::string, mycamera>& camera_cluster) {
//	std::vector<pts_on_mesh*> pts_dst;
//	for (int i = 0; i < pts_2d.size(); i++) {
//		ASSERT_WITH_MSG(pts_2d[i].size() == camera_cluster.size(), "The size of the input points and camera is not equal.");
//
//		pts_dst.push_back(pts_back_projection_multiview(pts_2d[i], camera_cluster));
//	}
//	return pts_dst;
//}

// deprecated
// M_mat is the camera K * Rt, debug mode
//pts_on_mesh* MyMesh::project_prediction(double x, double y, cv::Mat& M_mat, double conf, std::vector<double>& ray_found, std::vector<double>& camera_center, bool debug) {
//	std::vector<std::vector<double>> M(3, std::vector<double>(4, 0));
//	ASSERT_WITH_MSG(M_mat.type() == CV_64F, "The type of input projection matrix should be CV_64F.");
//	cv::Mat A(3, 3, CV_64F);
//	cv::Mat b(3, 1, CV_64F);
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 4; j++) {
//			M[i][j] = M_mat.at<double>(i, j);
//			if (j < 3) {
//				A.at<double>(i, j) = M[i][j];
//			}
//			else {
//				b.at<double>(i, 0) = M[i][3];
//			}
//		}
//	}
//
//	std::cout << "A is:" << std::endl;
//	print_mat(A);
//	std::cout << "B is:" << std::endl;
//	print_mat(b);
//	std::cout << "M_mat is:" << std::endl;
//	print_mat(M_mat);
//	cv::Mat center = -A.inv() * b;		// camera center coordinate
//
//	std::vector<double> C;
//	for (int i = 0; i < 3; i++) {
//		C.push_back(center.at<double>(i, 0));
//	}
//
//	// expand projection equation to find two planes
//	for (int i = 0; i < 3; i++) {
//		M[0][i] -= x * M[2][i];
//		M[1][i] -= y * M[2][i];
//	}
//	std::cout << "M is:" << std::endl;
//	for (int i = 0; i < M.size() - 1; i++) {
//		print_vec(M[i]);
//	}
//
//	// extract the first three variable from plane vector as the direction vector
//	std::vector<double> plane1, plane2;
//	for (int i = 0; i < 3; i++) {
//		plane1.push_back(M[0][i]);
//		plane2.push_back(M[1][i]);
//	}
//	std::vector<double> ray = cross(plane1, plane2);
//	std::cout << "initial ray is:" << std::endl;
//	print_vec(ray);
//	std::cout << std::endl;
//
//	double acc = 0;
//	for (int i = 0; i < 3; i++) {
//		acc += ray[i] * ray[i];
//	}
//	acc = sqrt(acc);
//	cv::Mat test_point(4, 1, CV_64F);
//	for (int i = 0; i < 3; i++) {
//		test_point.at<double>(i, 0) = C[i] + ray[i];
//	}
//	test_point.at<double>(3, 0) = 1;
//	test_point = M_mat * test_point;
//	if (test_point.at<double>(2, 0) < 0) {		//behind camera
//		acc *= -1;
//	}
//
//	std::cout << "dividor is:" << std::endl;
//	std::cout << acc << std::endl;
//	std::cout << std::endl;
//	for (int i = 0; i < 3; i++) {
//		ray[i] /= acc;
//	}
//
//
//	// output ray
//	for (int i = 0; i < 4; i++) {
//		ray_found.push_back(ray[i]);
//	}
//	// output camera center
//	for (int i = 0; i < 3; i++) {
//		camera_center.push_back(C[i]);
//	}
//	camera_center.push_back(1.0);
//	return this->get_pts_on_mesh(C, ray, conf);
//}

// deprecated
// M_mat is the camera K * Rt
//pts_on_mesh* MyMesh::project_prediction(double x, double y, cv::Mat& M_mat, double conf, std::vector<double>& ray_found, std::vector<double>& camera_center) {
//	std::vector<std::vector<double>> M(3, std::vector<double>(4, 0));
//	ASSERT_WITH_MSG(M_mat.type() == CV_64F, "The type of input projection matrix should be CV_64F.");
//	cv::Mat A(3, 3, CV_64F);
//	cv::Mat b(3, 1, CV_64F);
//	for (int i = 0; i < 3; i++) {
//		for (int j = 0; j < 4; j++) {
//			M[i][j] = M_mat.at<double>(i, j);
//			if (j < 3) {
//				A.at<double>(i, j) = M[i][j];
//			}
//			else {
//				b.at<double>(i, 0) = M[i][3];
//			}
//		}
//	}
//
//	cv::Mat center = -A.inv() * b;		// camera center coordinate
//
//	std::vector<double> C;
//	for (int i = 0; i < 3; i++) {
//		C.push_back(center.at<double>(i, 0));
//	}
//
//	// expand projection equation to find two planes
//	for (int i = 0; i < 3; i++) {
//		M[0][i] -= x * M[2][i];
//		M[1][i] -= y * M[2][i];
//	}
//
//	// extract the first three variable from plane vector as the direction vector
//	std::vector<double> plane1, plane2;
//	for (int i = 0; i < 3; i++) {
//		plane1.push_back(M[0][i]);
//		plane2.push_back(M[1][i]);
//	}
//	std::vector<double> ray = cross(plane1, plane2);
//
//	double acc = 0;
//	for (int i = 0; i < 3; i++) {
//		acc += ray[i] * ray[i];
//	}
//	acc = sqrt(acc);
//	cv::Mat test_point(4, 1, CV_64F);
//	for (int i = 0; i < 3; i++) {
//		test_point.at<double>(i, 0) = C[i] + ray[i];
//	}
//	test_point.at<double>(3, 0) = 1;
//	test_point = M_mat * test_point;
//	if (test_point.at<double>(2, 0) < 0) {		//behind camera
//		acc *= -1;
//	}
//
//	for (int i = 0; i < 3; i++) {
//		ray[i] /= acc;
//	}
//
//
//	// output ray
//	for (int i = 0; i < 4; i++) {
//		ray_found.push_back(ray[i]);
//	}
//	// output camera center
//	for (int i = 0; i < 3; i++) {
//		camera_center.push_back(C[i]);
//	}
//	camera_center.push_back(1.0);
//	return this->get_pts_on_mesh(C, ray, conf);
//}



// deprecated
// reproject 2d ray to 3d mesh and find the intersection point, debug mode
//pts_on_mesh* MyMesh::get_pts_on_mesh(std::vector<double>& C, std::vector<double>& ray, double conf, bool debug) {
//	double closest = -1;
//	std::vector<double> cpts;
//	int tri_id = -1; 
//	int planes = this->planes.size();
//	double t;
//	int k;
//	std::vector<double> pts_tmp;
//
//	std::cout << "final ray is:" << std::endl;
//	print_vec(ray);
//	std::cout << std::endl;
//
//	std::cout << "center is:" << std::endl;
//	print_vec(C);
//	std::cout << std::endl;
//
//	// go though all planes
//	for (k = 0; k < planes; k++) {
//		//if (!mesh.plane_projection_good[k]) {
//		//	continue;
//		//}
//
//		// plane is plane[0] * x + plane[1] * y + plane[2] * z + plane[3] = 0
//		// expected 3d point is [C[0] + t * ray[0], C[1] + t * ray[1], C[2] + t * ray[2]], t is the parameter we need to find
//		std::vector<double>& plane = this->planes[k];
//		double n = 0, d = 0;
//		for (int l = 0; l < 3; l++) {
//			n += C[l] * plane[l];
//			d += ray[l] * plane[l];
//		}
//		n += plane[3];
//		if (abs(d) < EPS) {
//			continue;
//		}
//		t = -n / d;
//		if (t < 0) {
//			continue;
//		}
//
//		pts_tmp.clear();
//		for (int l = 0; l < 3; l++) {
//			pts_tmp.push_back(t * ray[l] + C[l]);	// the 3d point
//		}
//		/*
//		std::cout << "plane is:" << std::endl;
//		print_vec(plane);
//		std::cout << std::endl;
//
//		std::cout << "k is:" << k << std::endl;
//		std::cout << std::endl;
//
//		std::cout << "t is:" << t << std::endl;
//		std::cout << std::endl;
//
//		std::cout << "n is:" << n << std::endl;
//		std::cout << std::endl;
//
//		std::cout << "d is:" << d << std::endl;
//		std::cout << std::endl;*/
//
//		std::vector<double> pts_append = pts_tmp;
//		pts_append.push_back(1.0);
//		//std::cout << "current point is:" << std::endl;
//		//print_vec(pts_append);
//		ASSERT_WITH_MSG(abs(inner(pts_append, plane)) < EPS, "Point found is not even on the extended triangle plane. The offset is " + std::to_string(abs(inner(pts_append, plane))));		// check the point with parameter t is inside the current plane
//
//		//std::cout << "id of 3 points around the plane are" << std::endl;
//		std::vector<int> vec_id = this->plane_pts_idx[k];
//		std::vector<double> tri_a, tri_b, tri_c;
//		tri_a.push_back(this->cloud->points[vec_id[0]].x);
//		tri_a.push_back(this->cloud->points[vec_id[0]].y);
//		tri_a.push_back(this->cloud->points[vec_id[0]].z);
//		tri_b.push_back(this->cloud->points[vec_id[1]].x);
//		tri_b.push_back(this->cloud->points[vec_id[1]].y);
//		tri_b.push_back(this->cloud->points[vec_id[1]].z);
//		tri_c.push_back(this->cloud->points[vec_id[2]].x);
//		tri_c.push_back(this->cloud->points[vec_id[2]].y);
//		tri_c.push_back(this->cloud->points[vec_id[2]].z);
//
//		//if (p[0] >= 0 && p[1] >= 0 && (p[0] + p[1] <= 1)) {	// check if the point is inside the triangle
//		if (point_triangle_test_3d(pts_tmp, tri_a, tri_b, tri_c)) {
//			if (tri_id == -1 || t < closest) {		// find the exact right plane which intersect with the ray and the intersection point is inside the plane.
//				closest = t;
//				tri_id = k;
//				cpts = pts_tmp;
//
//				//std::cout << "plane id is:" << tri_id << std::endl;
//				//std::cout << std::endl;
//				//std::cout << "plane is:" << std::endl;
//				//print_vec(this->planes[tri_id], 15);
//				//std::cout << std::endl;
//
//				//std::cout << "id of 3 points around the plane are" << std::endl;
//				//std::vector<int> ids = this->plane_pts_idx[tri_id];
//				//std::cout << std::endl;
//
//				//std::cout << "The 3d coordinate of that three points:" << std::endl;
//				//std::vector<double> pts1, pts2, pts3;
//				//pts1.push_back(this->cloud->points[ids[0]].x);
//				//pts1.push_back(this->cloud->points[ids[0]].y);
//				//pts1.push_back(this->cloud->points[ids[0]].z);
//				//pts2.push_back(this->cloud->points[ids[1]].x);
//				//pts2.push_back(this->cloud->points[ids[1]].y);
//				//pts2.push_back(this->cloud->points[ids[1]].z);
//				//pts3.push_back(this->cloud->points[ids[2]].x);
//				//pts3.push_back(this->cloud->points[ids[2]].y);
//				//pts3.push_back(this->cloud->points[ids[2]].z);
//				//print_vec(pts1, 15);
//				//print_vec(pts2, 15);
//				//print_vec(pts3, 15);
//				//std::cout << std::endl;
//
//				//std::cout << "temporary 3d points is" << std::endl;
//				//print_vec(cpts, 15);
//				//std::cout << std::endl;
//				//point_triangle_test_3d(cpts, pts1, pts2, pts3, 1);
//				//ASSERT_WITH_MSG(point_triangle_test_3d(cpts, pts1, pts2, pts3) == true, "Wrong!!!");
//			}
//		}
//	}
//
//
//	std::cout << "plane id is:" << tri_id << std::endl;
//	std::cout << std::endl;
//	std::cout << "plane is:" << std::endl;
//	print_vec(this->planes[tri_id]);
//	std::cout << std::endl;
//
//	std::cout << "id of 3 points around the plane are";
//	std::vector<int> ids = this->plane_pts_idx[tri_id];
//	std::cout << ids[0] << ", " << ids[1] << ", " << ids[2] << std::endl;
//	std::cout << std::endl;
//
//	std::cout << "The 3d coordinate of that three points:" << std::endl;
//	std::vector<double> pts1, pts2, pts3;
//	pts1.push_back(this->cloud->points[ids[0]].x);
//	pts1.push_back(this->cloud->points[ids[0]].y);
//	pts1.push_back(this->cloud->points[ids[0]].z);
//	pts2.push_back(this->cloud->points[ids[1]].x);
//	pts2.push_back(this->cloud->points[ids[1]].y);
//	pts2.push_back(this->cloud->points[ids[1]].z);
//	pts3.push_back(this->cloud->points[ids[2]].x);
//	pts3.push_back(this->cloud->points[ids[2]].y);
//	pts3.push_back(this->cloud->points[ids[2]].z);
//	print_vec(pts1);
//	print_vec(pts2);
//	print_vec(pts3);
//	std::cout << std::endl;
//
//	std::cout << "t is:" << t << std::endl;
//	std::cout << std::endl;
//	std::cout << "final 3d points is" << std::endl;
//	print_vec(cpts);
//	std::cout << std::endl;
//	//std::cout << "final 3d points is" << std::endl;
//	//print_vec(pts_tmp);
//	//std::cout << std::endl;
//
//	pts_on_mesh* ptr_mesh;
//	if (tri_id == -1) {
//		conf = 0;
//		fprintf(stderr, "no projected triangles!\n");
//
//		ptr_mesh = new pts_on_mesh(-1, 0.0, 0.0, 0.0, conf);
//		return ptr_mesh;
//	}
//
//	//ptr_mesh = new pts_on_mesh(tri_id, cpts[0], cpts[1], cpts[2], conf);
//	//ptr_mesh = new pts_on_mesh(this->plane_pts_idx[tri_id][0], cpts[0], cpts[1], cpts[2], conf);	// TODO: now it put a fake id and actually this point is not on the mesh, need another function to find it
//	//return ptr_mesh;
//	//cpts.push_back(conf);
//	pts_3d_conf pts_3d(cpts[0], cpts[1], cpts[2], conf);
//	return find_closest_pts_on_mesh(pts_3d, tri_id);
//}