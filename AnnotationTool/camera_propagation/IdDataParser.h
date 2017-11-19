#ifndef _ID_DATA_PARSER_H_
#define _ID_DATA_PARSER_H_

#include "myheader.h"

//#include "CameraModel2.h"
//#include "debug_print.h"

class mycamera;
class ColorCameraModel;

// the following functions should be compatible with Mugsy and Sociopticon

__declspec(dllexport) bool LoadIdCalibration(const std::string& calib_fpath, std::vector<mycamera>& out_cameras, const bool consider_skew = false, const bool consider_dist = true);

bool LoadIdColorImages(const std::string& color_image_fpath_template, const std::vector<ColorCameraModel>& color_cameras, const int frame_number, std::vector<cv::Mat>& out_color_images);

bool LoadIdDepthMaps(const std::string& depth_map_fpath_template, const std::vector<ColorCameraModel>& color_cameras, const int frame_number, std::vector<cv::Mat_<float>>& out_depth_maps);


// deprecated
//bool LoadIdCalibration(const std::string& calib_fpath, std::map<std::string, mycamera>& out_cameras, const bool consider_skew = false, const bool consider_dist = true);


#endif  // _ID_DATA_PARSER_H_