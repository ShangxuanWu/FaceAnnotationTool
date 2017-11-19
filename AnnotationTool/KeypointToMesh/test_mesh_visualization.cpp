// Author: Xinshuo
// Email: xinshuow@andrew.cmu.edu


#include "mesh_visualization.h"


// pcl library
#include <pcl/io/auto_io.h>

int main() {
	// define mesh and get cloud
	pcl::PolygonMesh::Ptr mesh(new pcl::PolygonMesh);
	const std::string filename = "simple.ply";
	pcl::io::load(filename, *mesh.get());
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	fromPCLPointCloud2(mesh->cloud, *cloud.get());

	// define arbitraty points
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr keypoints_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	std::vector<std::vector<double>> points;
	std::vector<double> single_point;
	single_point.push_back(161.68);
	single_point.push_back(288.36);
	single_point.push_back(34.37);
	points.push_back(single_point);
	get_cloud_from_points(points, keypoints_cloud_ptr);

	//// define a arbitrary line
	//std::vector<double> line;
	//line.push_back(1);
	//line.push_back(1);
	//line.push_back(1);
	//line.push_back(1);
	//pcl::PointCloud<pcl::PointXYZRGB>::Ptr line_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);
	//get_cloud_from_line(line, line_cloud_ptr);

	// define a set of arbitrary line
	std::vector<std::vector<float>> lines;
	std::vector<float> line;
	line.push_back(1);
	line.push_back(1);
	line.push_back(1);
	line.push_back(0);
	lines.push_back(line);
	line.clear();
	line.push_back(1);
	line.push_back(1);
	line.push_back(-1);
	line.push_back(0);
	lines.push_back(line);
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr line_cloud_ptr(new pcl::PointCloud<pcl::PointXYZRGB>);

	std::vector<pcl::PointXYZ> pts_start;
	pts_start.push_back(pcl::PointXYZ(0, 0, 0));
	pts_start.push_back(pcl::PointXYZ(0, 0, 0));
	std::vector<uint32_t> range;
	range.push_back(64);
	range.push_back(128);
	get_cloud_from_lines(lines, line_cloud_ptr, pts_start, range);

	// test function
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer;
	//viewer = keypoint_mesh_visualization(keypoints_cloud_ptr, cloud);
	//std::map<std::string, pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr> cloud_map;
	//cloud_map["mesh"] = cloud;
	//cloud_map["keypoints"] = keypoints_cloud_ptr;
	//cloud_map["line"] = line_cloud_ptr;
	//viewer = rgb_cloud_visualization(cloud_map);
	viewer = keypoint_line_mesh_visualization(keypoints_cloud_ptr, line_cloud_ptr, cloud);
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}

}

