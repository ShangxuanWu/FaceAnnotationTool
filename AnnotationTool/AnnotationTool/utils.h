#ifndef UTILS_H
#define UTILS_H

#define PI 3.14159265

#define ORIGINAL_IMAGE_WIDTH 5120
#define ORIGINAL_IMAGE_HEIGHT 3840

// for image rotation
#define MIDDLE 0
#define LEFT -1
#define RIGHT 1

#define TOTAL_VIEW_COUNT 40
#define TOTAL_POINTS_IN_SINGLE_VIEW 68
#define DLIB_MODEL_PATH "shape_predictor_68_face_landmarks.dat"

#define QPAINTER_FONT_SIZE 55
#define MAX_SELECT_DISTANCE 8000
#define MAX_SELECT_DISTANCE_3D 10
#define PREPROCESS_RANSAC_ITER 2000
#define PROPAGATE_RANSAC_ITER 100

#define INITIAL_CONF 0.5
#define CONF_AFTER_FIRST_PROPAGATION_IN_PREPROCESS 0.0

// for the render mode in 3D panel
#define POINT_CLOUD_MODE 1
#define MESH_MODE 2

#define NOW_RENDER_MODE MESH_MODE
#define POINT_SIZE_3D 50
#define POINT_MINIMAL_SIZE_3D 1

// for the instruction panel
#define INSTRUCTION_IMG_WIDTH 603
#define INSTRUCTION_IMG_HEIGHT 442

// for the Epipolar Window
#define EPIPOLAR_HEIGHT 600 
#define EPIPOLAR_WIDTH 600

// for the two small windows in epipolar window
#define EPIPOLAR_LEFT_WINDOW 1
#define EPIPOLAR_RIGHT_WINDOW 2

#define FIND_POINT_VISIBILE_TOLERANCE 1

#include <opencv2/opencv.hpp>

#include <unordered_map>
#include <string>
#include <vector>
#include <utility>

#include <QDebug>
#include <QString>
#include <set>

#include "../camera_propagation/pts_2d_tool.h"
#include "../camera_propagation/mycamera.h"

// which direction the image <want> to rotate in order to get an upright image
static std::unordered_map<std::string, int> rotation_indicator{
{"01",LEFT },
{"05",RIGHT},
{"06",RIGHT},
{"07",RIGHT},
{"10",LEFT},
{"11",RIGHT},
{"12",RIGHT},
{"13",RIGHT},
{"14",RIGHT},
{"15",RIGHT},
{"16",RIGHT},
{"17",RIGHT},
{"18",RIGHT},
{"19",LEFT},
{"20",LEFT},
{"21",LEFT},
{"22",LEFT},
{"23",RIGHT},
{"24",RIGHT},
{"25",RIGHT},
{"26",RIGHT},
{"27",RIGHT},
{"28",RIGHT},
{"29",RIGHT},
{"30",RIGHT},
{"31",RIGHT},
{"32",RIGHT},
{"33",RIGHT},
{"34",RIGHT},
{"35",LEFT},
{"36",RIGHT},
{"37",RIGHT},
{"38",LEFT},
{"39",RIGHT},
{"40",RIGHT},
{"41",RIGHT},
{"42",RIGHT},
{"43",RIGHT},
{"44",RIGHT},
{"45",RIGHT},
};

// for drawing connection lines in 68 points
static std::unordered_map<int, int> connected_points{
	// cheek
	{1,2},
	{ 2,3 },
	{ 3,4 },
	{ 4,5 },
	{ 5,6 },
	{ 6,7 },
	{ 7,8 },
	{ 8,9 },
	{ 9,10 },
	{ 10,11 },
	{ 11,12 },
	{ 12,13 },
	{ 13,14 },
	{ 14,15 },
	{ 15,16 },
	{ 16,17 },
	// left eyebrow
	{ 18,19 },
	{ 19,20 },
	{ 20,21 },
	{ 21,22 },
	// left eye
	{ 37,38 },
	{ 38,39 },
	{ 39,40 },
	{ 40,41 },
	{ 41,42 },
	{ 42,37 },
	// right eyebrow
	{ 23,24 },
	{ 24,25 },
	{ 25,26 },
	{ 26,27 },
	// right eye
	{ 43,44 },
	{ 44,45 },
	{ 45,46 },
	{ 46,47 },
	{ 47,48 },
	{ 48,43 },
	// nose
	{ 28,29 },
	{ 29,30 },
	{ 30,31 },
	// nose bottom part
	{ 32,33 },
	{ 33,34 },
	{ 34,35 },
	{ 35,36 },
	// mouth outer
	{ 49,50 },
	{ 50,51 },
	{ 51,52 },
	{ 52,53 },
	{ 53,54 },
	{ 54,55 },
	{ 55,56 },
	{ 56,57 },
	{ 57,58 },
	{ 58,59 },
	{ 59,60 },
	{ 60,49 },
	// mouth inner
	{ 61,62 },
	{ 62,63 },
	{ 63,64 },
	{ 64,65 },
	{ 65,66 },
	{ 66,67 },
	{ 67,68 },
	{ 68,61 },
};

static std::unordered_map<int, std::pair<std::string, std::string>> point_to_epipolar_image_choice{
	// left cheek
	{ 1,{ "330022","330007" } },
	{ 2,{ "330022","330007" } },
	{ 3,{ "330022","330007" } },
	{ 4,{ "330022","330007" } },
	{ 5,{ "330022","330007" } },
	{ 6,{ "330022","330007" } },
	{ 7,{ "330022","330007" } },

	// center cheek
	{ 8,{ "330010","330012" } },
	{ 9,{ "330010","330012" } },
	{ 10,{ "330010","330012" } },
	
	// right cheek
	{ 11,{ "330025","330035" } },
	{ 12,{ "330025","330035" } },
	{ 13,{ "330025","330035" } },
	{ 14,{ "330025","330035" } },
	{ 15,{ "330025","330035" } },
	{ 16,{ "330025","330035" } },
	{ 17,{ "330025","330035" } },
	
	// left eyebrow
	{ 18,{ "330016","330034" } },
	{ 19,{ "330016","330034" } },
	{ 20,{ "330016","330034" } },
	{ 21,{ "330016","330034" } },
	{ 22,{ "330016","330034" } },
		
	// right eyebrow
	{ 23,{ "330028","330045" } },
	{ 24,{ "330028","330045" } },
	{ 25,{ "330028","330045" } },
	{ 26,{ "330028","330045" } },
	{ 27,{ "330028","330045" } },
	
	// nose
	{ 28,{ "330016","330028" } },
	{ 29,{ "330016","330028" } },
	{ 30,{ "330016","330028" } },
	{ 31,{ "330016","330028" } },
	
	// bottom of nose
	{ 32,{ "330010","330012" } },
	{ 33,{ "330010","330012" } },
	{ 34,{ "330010","330012" } },
	{ 35,{ "330010","330012" } },
	{ 36,{ "330010","330012" } },
	
	// left eye
	{ 37,{ "330016","330034" } },
	{ 38,{ "330016","330034" } },
	{ 39,{ "330016","330034" } },
	{ 40,{ "330016","330034" } },
	{ 41,{ "330016","330034" } },
	{ 42,{ "330016","330034" } },
	
	// right eye
	{ 43,{ "330028","330045" } },
	{ 44,{ "330028","330045" } },
	{ 45,{ "330028","330045" } },
	{ 46,{ "330028","330045" } },
	{ 47,{ "330028","330045" } },
	{ 48,{ "330028","330045" } },

	// mouth
	{ 49,{ "330014","330031" } },
	{ 50,{ "330014","330031" } },
	{ 51,{ "330014","330031" } },
	{ 52,{ "330014","330031" } },
	{ 53,{ "330014","330031" } },
	{ 54,{ "330014","330031" } },
	{ 55,{ "330014","330031" } },
	{ 56,{ "330014","330031" } },
	{ 57,{ "330014","330031" } },
	{ 58,{ "330014","330031" } },
	{ 59,{ "330014","330031" } },
	{ 60,{ "330014","330031" } },
	{ 61,{ "330014","330031" } },
	{ 62,{ "330014","330031" } },
	{ 63,{ "330014","330031" } },
	{ 64,{ "330014","330031" } },
	{ 65,{ "330014","330031" } },
	{ 66,{ "330014","330031" } },
	{ 67,{ "330014","330031" } },
	{ 68,{ "330014","330031" } },
};

// for drawing connection lines in 192 points
//static std::unordered_map<int, int> connected_points{}

// this is for 40 views, the order of camera in the right panel list
static std::vector<std::string> camera_order = {"19", "05", "32", "17", "11", "16", "45", "23", "24", "36", "39", "22", "27", "44", "26", "18", "15", "06", "34", "14", "30", "28", "33", "43", "41", "37", "21", "40", "35", "07", "13", "29", "01", "10", "31", "12", "20", "38", "42", "25"};

static std::set<std::string> valid_dlib_view = { "330038", "330010", "330001", "330020" };

// get the last 2 char of the folder name as ID
std::string getIDFromFilename(std::string str);

// get 3300XX from ID
std::string getFolderFromFilename(std::string full_path);

std::string getIDFromFoldername(std::string folder_name);

std::string getFolderFromID(std::string ID);

std::string getFullPathFromFolder(std::string root, std::string folder_name);

// used in onPopUpEpipolarWindow()
mycamera getMycameraFromFullPath(QString imagePaths, std::vector<mycamera> camera_cluster);


void rot90(cv::Mat &matImage, int rotflag);

// helper for visualizing the points
void helper_map_to_vec_tool_to_conf(std::map<std::string, std::vector<pts_2d_conf>> pts_src);

#endif