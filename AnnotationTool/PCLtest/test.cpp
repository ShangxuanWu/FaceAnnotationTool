#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <iostream>
#include <opencv/cv.h>
#include <pcl/io/auto_io.h>
//#include <pcl/io/ply_io.h>
//#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/PolygonMesh.h>
#include <pcl/conversions.h>
#include <fstream>
#include <unordered_set>

#include <pcl/io/boost.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/io/ifs_io.h>
#include <pcl/io/obj_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/io/ply/ply_parser.h>


using namespace std;
using namespace cv;
using namespace pcl;
using namespace pcl::io;
#define EPS 1e-10



int main(int argc, char* argv[]) {
	//if (argc != 8) {
	//	fprintf(stderr, "number of arguments is %d", argc);
	//	fprintf(stderr, "calibration pose_dir mesh num_parts frame pose3dout_dir resize_factor\n");
	//	return -1;
	//}
	//char filename1[1000] = "D:\oculus\NewAnnotationTool\AnnotationTool_Oculus\KeypointToMesh\00000.ply";
	//char filename1[1000] = "D:\oculus\NewAnnotationTool\AnnotationTool_Oculus\KeypointToMesh\00000.ply";
	char filename[1000] = "simple.ply";
	//int scale = 1;

	pcl::PointCloud<pcl::PointXYZ> cloud;
	pcl::PolygonMesh mesh;
	//pcl::PolygonMesh mesh;
	const std::string file_temp = filename;
	//pcl::PolygonMesh mesh;
	//pcl::io::load(file_temp, cloud);
	//pcl::io::loadPolygonFileVTK("simple.ply", mesh);
	//pcl::io::load(file_temp, *mesh.get());
	//pcl::io::load(filename, *mesh.get());
	//pcl::io::loadPolygonFilePLY(filename1, *mesh);
	//pcl::io::loadPolygonFilePLY(file_temp, *mesh);
	//pcl::io::loadPLYFile(filename1, *mesh.get());
	//pcl::io::loadPLYFile(file_temp, *mesh.get());

	pcl::PCLPointCloud2 cloud_blob, cloud_blob2;
	PLYReader reader;
	reader.read(file_temp, cloud_blob2);



	//cloud = pcl::PointCloud<pcl::PointXYZ>::Ptr(new pcl::PointCloud<pcl::PointXYZ>);
	//fromPCLPointCloud2(mesh->cloud, *cloud.get());
	//for (unsigned int i = 0; i < cloud->size(); i++) {
	//	cloud->points[i].x *= scale;
	//	cloud->points[i].y *= scale;
	//	cloud->points[i].z *= scale;
	//}




	return 0;
}





