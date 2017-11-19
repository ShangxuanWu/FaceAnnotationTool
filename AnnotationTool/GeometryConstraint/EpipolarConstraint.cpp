#include "EpipolarConstraint.h"
//#include <opencv2/core/eigen.hpp>
#include <Eigen/LU>
#include <Eigen/SVD>
#include <Eigen/Core>

EpipolarConstraint::EpipolarConstraint(std::string root_folder_)
{
	root_folder = root_folder_;

	// read the camera parameters
	std::map<std::string, std::vector<pts_2d_conf>> pts_src;
	std::vector<mycamera> camera_cluster;
	std::map<std::string, std::vector<pts_2d_conf>> pts_output;

	// list the dirs and read the files
	std::string calib_txt_path = root_folder + "\\calibration.txt";
	getSubDirs(camera_names, root_folder);
	readCameraMatrix(calib_txt_path, cameras);

	// initialize nums and cameras
	num1 = 0;
	num2 = 1;
	camera1_id = findCameraByName(camera_names[num1].substr(3));
	camera2_id = findCameraByName(camera_names[num2].substr(3));

	clicked_x = -1;
	clicked_y = -1;

}

int EpipolarConstraint::findCameraByName(std::string camera_name)
{
	for (int i = 0; i < cameras.size(); i++)
	{
		if (camera_name == cameras[i].name)
		{
			//return cameras[i];
			return i;
		}
	}
}

void EpipolarConstraint::onMouse(int event, int x, int y, int, void* userdata)
{
	// Check for null pointer in userdata and handle the error
	EpipolarConstraint* epipolar_constraint = reinterpret_cast<EpipolarConstraint*>(userdata);
	epipolar_constraint->onMouse(event, x, y);
}

void EpipolarConstraint::onMouse(int event, int x, int y)
{
	if (event == cv::EVENT_LBUTTONDOWN)
	{
		std::cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << std::endl;
		// goto calculate epipolar line
		clicked_x = x * downsize_ratio;
		clicked_y = y * downsize_ratio;
	}
}

void EpipolarConstraint::getClickedPoint()
{
	std::string img1_fn = root_folder + '\\' + camera_names[num1] + "\\image0000.png";
	std::string img2_fn = root_folder + '\\' + camera_names[num2] + "\\image0000.png";

	
	image1 = cv::imread(img1_fn, CV_LOAD_IMAGE_COLOR);
	image2 = cv::imread(img2_fn, CV_LOAD_IMAGE_COLOR);

	cv::namedWindow("Image1", cv::WINDOW_AUTOSIZE);// Create a window for display.
	cv::setMouseCallback("Image1", onMouse, this);

	size = cv::Size(int(image1.cols / downsize_ratio), int(image1.rows / downsize_ratio));

	cv::resize(image1, image1_resized, size);
	//cv::resize(image2, image2_resized, size);


	cv::imshow("Image1", image1_resized);                   // Show our image inside it.
	//cv::namedWindow("Image2", cv::WINDOW_AUTOSIZE);// Create a window for display.
	//cv::imshow("Image2", image2_resized);                   // Show our image inside it.

	cv::waitKey(0);

}

void EpipolarConstraint::calculateEpipolarLine()
{
	std::cout << "Calculating Epipolar Line !" << std::endl;
	calculateEpipolarLine(cv::Point(clicked_x, clicked_y));
}

void EpipolarConstraint::showEpipolarLine()
{
	cv::Size size2 = cv::Size(int(image1.cols / downsize_ratio2), int(image1.rows / downsize_ratio2));
	cv::resize(image2, image2_resized, size2);
	
	std::cout << "---line---" << std::endl;
	std::cout << line[0] << std::endl;
	std::cout << line[1] << std::endl;
	std::cout << line[2] << std::endl;
	std::cout << "---points---" << std::endl;
	for (int x = 0; x < image2.cols; x++)
	{
		int y = -(line[2] + line[0]*x) / line[1];
		if (y >= 0 && y < image2.rows)
		{
			int x1 = x / downsize_ratio2;
			int y1 = y / downsize_ratio2;
			std::cout << x << ' ' << y << std::endl;
			//image2.at<cv::Vec3b>(y, x) = 255;
			image2_resized.at<cv::Vec3b>(y1, x1)[0] = 255;
			image2_resized.at<cv::Vec3b>(y1, x1)[1] = 255;
			image2_resized.at<cv::Vec3b>(y1, x1)[2] = 255;
			//image2.at<uchar>(y, x, 0) = 255;
			//image2.at<uchar>(y, x, 1) = 255;
			//image2.at<uchar>(y, x, 2) = 255;
			// draw white line
		}
	}
	//cv::imwrite("result.png", image2);
	cv::namedWindow("Image2", cv::WINDOW_AUTOSIZE);// Create a window for display.
	cv::imshow("Image2", image2_resized);                   // Show our image inside it.
	//cv::imshow("Image2", image2);                   // Show our image inside it.
	cv::waitKey(0);
}

void EpipolarConstraint::readCameraMatrix(std::string calibration_file_path, std::vector<mycamera>& camera_names)
{
	LoadIdCalibration(calibration_file_path, cameras, false);
}

/// Gets a list of subdirectories under a specified path
/// @param[out] output Empty vector to be filled with result
/// @param[in]  path   Input path, may be a relative path from working dir
void EpipolarConstraint::getSubDirs(std::vector<std::string>& output, const std::string& path)
{
	// clear output first
	output.clear();

	WIN32_FIND_DATA findfiledata;
	HANDLE hFind = INVALID_HANDLE_VALUE;

	char fullpath[MAX_PATH];
	GetFullPathName(path.c_str(), MAX_PATH, fullpath, 0);
	std::string fp(fullpath);

	hFind = FindFirstFile((LPCSTR)(fp + "\\*").c_str(), &findfiledata);
	if (hFind != INVALID_HANDLE_VALUE)
	{
		do
		{
			if ((findfiledata.dwFileAttributes | FILE_ATTRIBUTE_DIRECTORY) == FILE_ATTRIBUTE_DIRECTORY
				&& (findfiledata.cFileName[0] != '.'))
			{
				output.push_back(findfiledata.cFileName);
			}
		} while (FindNextFile(hFind, &findfiledata) != 0);
	}
}

cv::Mat getEpipolarLine(mycamera camera1, mycamera camera2, cv::Point original_point1)
{
	cv::Mat point1 = (cv::Mat_<double>(3, 1) << original_point1.x, original_point1.y, 1);
	cv::Mat P1 = camera1.intrinsic * camera1.extrinsic;
	cv::Mat P2 = camera2.intrinsic * camera2.extrinsic;
	// how to get null space? ///////////////////////////
	cv::Mat C = getNullSpace(P1);
	print_cv_mat(C);
	///////////////////////////////////////////////////////
	//std::cout << C.type() << std::endl;
	//std::cout << P2.type() << std::endl;
	cv::Mat e_ = P2 * C;
	cv::Mat X = (cv::Mat_<double>(3, 3) << 0, -e_.at<double>(2, 0), e_.at<double>(1, 0), e_.at<double>(2, 0), 0, -e_.at<double>(0, 0), -e_.at<double>(1, 0), e_.at<double>(0, 0), 0);
	//std::cout << X << std::endl;
	cv::Mat line_ = X * P2 * P1.inv(cv::DECOMP_SVD) * point1;
	return line_;
}

void EpipolarConstraint::calculateEpipolarLine(cv::Point original_point1)
{
	if (clicked_x == -1 || clicked_y == -1)
	{
		std::cout << "Didn't get the point! Exiting!" << std::endl;
	}
	
	// call external function getEpipolarLine()
	cv::Mat line_ = getEpipolarLine(cameras[camera1_id], cameras[camera2_id], original_point1);
	/////////////////////////////////////////////

	//std::cout << line_.rows;
	//std::cout << line_.cols;
	//return line_;
	line[0] = line_.at<double>(0, 0);
	line[1] = line_.at<double>(1, 0);
	line[2] = line_.at<double>(2, 0);
}

// use Eigen library to get null space here
cv::Mat getNullSpace(cv::Mat p)
{
	cv::SVD svd = cv::SVD(p, cv::SVD::FULL_UV);
	cv::Mat vt_ = svd.vt;
	cv::Mat u = svd.u;
	cv::Mat w = svd.w;

	int i;
	for (i = 1; i <= 3; i++)
	{
		if (p.at<double>(i - 1, i - 1) == 0)
		{
			break;
		}
	}

	cv::Mat result = vt_(cv::Rect(0, i-1, p.cols, vt_.rows-i+1));
	cv::Mat result_t;
	cv::transpose(result, result_t);
	return result_t;
}

void print_cv_mat(cv::Mat P)
{
	std::cout << "---Printing Matrix---" << std::endl;
	for (int i = 0; i < P.rows; i++)
	{
		for (int j = 0; j < P.cols; j++)
		{
			std::cout << P.at<double>(i, j) << ' ';
		}
		std::cout << std::endl;
	}
	std::cout << "---------------------" << std::endl;
}