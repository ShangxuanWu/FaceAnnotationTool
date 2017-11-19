#ifndef MAIN_FUNCTION_H
#define MAIN_FUNCTION_H

#define ENABLE_ASSERTS

#include <opencv2/opencv.hpp>
#include <opencv2/core/core_c.h>
#include <iostream>
#include <string>
#include <vector>
#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/image_processing/render_face_detections.h>
#include <dlib/image_processing.h>
#include <dlib/opencv/cv_image.h>
#include <dlib/gui_widgets.h>
#include <dlib/image_io.h>

__declspec(dllexport) std::vector<cv::Point> detect_one_image_using_dlib(std::string file_path, std::string model_path);
__declspec(dllexport) std::vector<cv::Point> detect_one_image_using_dlib_input_Mat(const cv::Mat& img, std::string model_path);

#endif // MAIN_FUNCTION_H