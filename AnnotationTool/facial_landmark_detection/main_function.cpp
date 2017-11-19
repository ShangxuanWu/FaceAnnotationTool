// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

This example program shows how to find frontal human faces in an image and
estimate their pose.  The pose takes the form of 68 landmarks.  These are
points on the face such as the corners of the mouth, along the eyebrows, on
the eyes, and so forth.



This face detector is made using the classic Histogram of Oriented
Gradients (HOG) feature combined with a linear classifier, an image pyramid,
and sliding window detection scheme.  The pose estimator was created by
using dlib's implementation of the paper:
One Millisecond Face Alignment with an Ensemble of Regression Trees by
Vahid Kazemi and Josephine Sullivan, CVPR 2014
and was trained on the iBUG 300-W face landmark dataset.

Also, note that you can train your own models using dlib's machine learning
tools.  See train_shape_predictor_ex.cpp to see an example.




Finally, note that the face detector is fastest when compiled with at least
SSE2 instructions enabled.  So if you are using a PC with an Intel or AMD
chip then you should enable at least SSE2 instructions.  If you are using
cmake to compile this program you can enable them by using one of the
following commands when you create the build project:
cmake path_to_dlib_root/examples -DUSE_SSE2_INSTRUCTIONS=ON
cmake path_to_dlib_root/examples -DUSE_SSE4_INSTRUCTIONS=ON
cmake path_to_dlib_root/examples -DUSE_AVX_INSTRUCTIONS=ON
This will set the appropriate compiler options for GCC, clang, Visual
Studio, or the Intel compiler.  If you are using another compiler then you
need to consult your compiler's manual to determine how to enable these
instructions.  Note that AVX is the fastest but requires a CPU from at least
2011.  SSE4 is the next fastest and is supported by most current machines.
*/

#define ENABLE_ASSERTS

#include "main_function.h"

// ----------------------------------------------------------------------------------------

std::vector<cv::Point> detect_one_image_using_dlib(std::string file_path, std::string model_path)
{
	cv::Mat img = cv::imread(file_path, CV_LOAD_IMAGE_COLOR);
	dlib::array2d<dlib::rgb_pixel> img_dlib_format;
	dlib::assign_image(img_dlib_format, dlib::cv_image<dlib::bgr_pixel>(img));

	std::vector<cv::Point> result;
	dlib::shape_predictor sp;
	dlib::deserialize(model_path) >> sp;
	//dlib::rectangle rec(10,10, img.cols - 10, img.rows - 10);
	dlib::frontal_face_detector detector = dlib::get_frontal_face_detector();
	std::vector<dlib::rectangle> dets = detector(img_dlib_format);
	// if we do not use detector first, the program will be a lot faster
	
	if (dets.size() > 0) {
		dlib::rectangle rec = dets[0];
		dlib::full_object_detection shape = sp(img_dlib_format, rec);
		for (int k = 0; k < shape.num_parts(); k++) {
			cv::Point landmark(shape.part(k).x(),
				shape.part(k).y());
			result.push_back(landmark);
		}
	}
	return result;
}

// ----------------------------------------------------------------------------------------

std::vector<cv::Point> detect_one_image_using_dlib_input_Mat(const cv::Mat& img, std::string model_path)
{
	//cv::Mat img = cv::imread(file_path, CV_LOAD_IMAGE_COLOR);
	dlib::array2d<dlib::rgb_pixel> img_dlib_format;
	dlib::assign_image(img_dlib_format, dlib::cv_image<dlib::bgr_pixel>(img));

	std::vector<cv::Point> result;
	dlib::shape_predictor sp;
	dlib::deserialize(model_path) >> sp;
	dlib::frontal_face_detector detector = dlib::get_frontal_face_detector();
	std::vector<dlib::rectangle> dets = detector(img_dlib_format);

	if (dets.size() > 0) {
		dlib::rectangle rec = dets[0];
		dlib::full_object_detection shape = sp(img_dlib_format, rec);
		for (int k = 0; k < shape.num_parts(); k++) {
			cv::Point landmark(shape.part(k).x(),
				shape.part(k).y());
			result.push_back(landmark);
		}
	}
	return result;
}