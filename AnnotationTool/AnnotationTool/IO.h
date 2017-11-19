#ifndef IO_H
#define IO_H

// for saving and loading xml files storing all the annotated points

#include <iostream>
#include <string>
#include <map>

#include "opencv2/opencv.hpp"

#include "utils.h"
#include "../camera_propagation/pts_2d_tool.h"


void save_annotated_points(std::string filename, std::map<std::string, std::map<int, pts_2d_tool>>);
void save_points_visibility(std::string filename, std::map<std::string, std::vector<bool>> total_view_visibility);

std::map<std::string, std::map<int, pts_2d_tool>> load_annotated_points(std::string filename);

#endif