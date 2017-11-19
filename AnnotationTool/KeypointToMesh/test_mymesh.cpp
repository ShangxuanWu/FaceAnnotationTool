// Modified by: Xinshuo
// Email: xinshuow@andrew.cmu.edu

//#include "../camera_propagation/io_point.h"
#include "../camera_propagation/debug_print.h"

// in-project library
#include "myheader.h"
#include "MyMesh.h"
//#include "special_camera.h"
//#include "3d_geometry.h"
//#include "math_functions.h"
//#include "pts_on_mesh.h"
#include "mesh_visualization.h"

int main(int argc, char* argv[]) {
	std::cout << "reading mesh file" << std::endl;
	const clock_t begin_time = clock();
	MyMesh mesh("sequence\\input\\meshes\\00000.ply", 1);			// input reconstructed mesh
	std::cout << "It spend " << float(clock() - begin_time) / CLOCKS_PER_SEC << " second to read the mesh file" << std::endl;

	int plane_id = 3473258;
	// query points id from plane id 
	std::cout << "Id of points corresponds to the first plane:" << std::endl;
	std::vector<int> ids = mesh.plane_pts_idx[plane_id];
	ASSERT_WITH_MSG(ids.size() == 3, "The size of id corresponds to first plane is not correct!");
	print_vec(ids);

	// print points coordinate
	std::cout << "The 3d coordinate of that three points:" << std::endl;
	std::vector<double> pts1, pts2, pts3;
	pts1.push_back(mesh.cloud->points[ids[0]].x);
	pts1.push_back(mesh.cloud->points[ids[0]].y);
	pts1.push_back(mesh.cloud->points[ids[0]].z);
	pts2.push_back(mesh.cloud->points[ids[1]].x);
	pts2.push_back(mesh.cloud->points[ids[1]].y);
	pts2.push_back(mesh.cloud->points[ids[1]].z);
	pts3.push_back(mesh.cloud->points[ids[2]].x);
	pts3.push_back(mesh.cloud->points[ids[2]].y);
	pts3.push_back(mesh.cloud->points[ids[2]].z);
	print_vec(pts1);
	print_vec(pts2);
	print_vec(pts3);
	
	// query plane id from points id
	std::vector<int> plane_id1 = mesh.planes_of_vertices[ids[0]];
	std::vector<int> plane_id2 = mesh.planes_of_vertices[ids[1]];
	std::vector<int> plane_id3 = mesh.planes_of_vertices[ids[2]];
	std::vector<int>::iterator it = std::find(plane_id1.begin(), plane_id1.end(), plane_id);
	std::cout << "Plane the first point belongs to is" << std::endl;
	print_vec(plane_id1);
	std::cout << "Plane the second point belongs to is" << std::endl;
	print_vec(plane_id2); 
	std::cout << "Plane the third point belongs to is" << std::endl;
	print_vec(plane_id3);
	ASSERT_WITH_MSG(it < plane_id1.end(), "Plane not found.");
	it = std::find(plane_id2.begin(), plane_id2.end(), plane_id);
	ASSERT_WITH_MSG(it < plane_id2.end(), "Plane not found.");
	it = std::find(plane_id3.begin(), plane_id3.end(), plane_id);
	ASSERT_WITH_MSG(it < plane_id3.end(), "Plane not found.");
	
	system("PAUSE");

	return 0;
}