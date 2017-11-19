#include "IO.h"
#include <QDebug>

void save_annotated_points(std::string filename, std::map<std::string, std::map<int, pts_2d_tool>> all_annotated_points) 
{
	cv::FileStorage fs(filename, cv::FileStorage::WRITE);

	// save camera names
	fs << "Cameras" << "[";

	for (auto it = all_annotated_points.begin(); it != all_annotated_points.end(); it++) {
		fs << "{:" << "FolderName" << it->first << "}";
	}
	fs << "]";
	//
	
    // saving points
	fs << "Points" << "[";
	
	for (auto it = all_annotated_points.begin(); it != all_annotated_points.end(); it++) {
		for (auto it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			fs << "{:" <<"CameraLabel" << it->first << "PointLabel" << it2->first << "X" << it2->second.x << "Y" << it2->second.y << "Anchor" << it2->second.anchor << "Confidence" << it2->second.conf << "}";
		}
	}
	fs << "]";
	//

	fs.release();
}

void save_points_visibility(std::string filename, std::map<std::string, std::vector<bool>> total_view_visibility)
{
	cv::FileStorage fs(filename, cv::FileStorage::WRITE);

	// save camera names
	fs << "Cameras" << "[";

	for (auto it = total_view_visibility.begin(); it != total_view_visibility.end(); it++) {
		fs << "{:" << "FolderName" << it->first << "}";
	}
	fs << "]";
	//

	// saving point visibilities
	fs << "Points" << "[";

	for (auto it = total_view_visibility.begin(); it != total_view_visibility.end(); it++) {
		for (int it2 = 0; it2 < TOTAL_POINTS_IN_SINGLE_VIEW; it2++) {
			fs << "{:" << "CameraLabel" << it->first << "PointLabel" << it2+1 << "Visibility" << it->second[it2] << "}";
		}
	}
	fs << "]";
	//

	fs.release();
}

std::map<std::string, std::map<int, pts_2d_tool>> load_annotated_points(std::string filename) {
	cv::FileStorage fs(filename, cv::FileStorage::READ);

	std::map<std::string, std::map<int, pts_2d_tool>> result;
	
	// retrieving camera names
	// iterate through a sequence using FileNodeIterator
	cv::FileNode cameras = fs["Cameras"];
	cv::FileNodeIterator it = cameras.begin(), it_end = cameras.end();
	for (; it != it_end; ++it)
	{
		std::string folder_name = (std::string)(*it)["FolderName"];
		std::map<int, pts_2d_tool> now_map;
		result[folder_name] = now_map;
	}

	// retrieving points
	// iterate through a sequence using FileNodeIterator
	cv::FileNode points = fs["Points"];
	it = points.begin(), it_end = points.end();
	for (; it != it_end; ++it)
	{
		std::string camera_label;
		int point_label;
		float	x, y, conf;
		bool anchor;
		camera_label = (std::string)(*it)["CameraLabel"];
		point_label = (int)(*it)["PointLabel"];
		x = (float)(*it)["X"];
		//qDebug() << x;
		y = (float)(*it)["Y"];
		//qDebug() << y;
		conf = (float)(*it)["Confidence"];
		//qDebug() << conf;
		anchor = (int)(*it)["Anchor"];
		
		result[camera_label][point_label] = pts_2d_tool(x, y, conf, anchor);
	}
	fs.release();
	return result;
}