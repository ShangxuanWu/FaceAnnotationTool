// Author: Xinshuo Weng
// Email: xinshuow@andrew.cmu.edu
#include "CameraModel2.h"

#include <iostream>
#include <fstream>

//#include "CalibrationIO.h"
//#include "DepthMapIO.h"

//#include "ListFileIO.h"

using namespace std;

CameraModel::CameraModel()
{
	globalOrient_mat_ = Eigen::Matrix3d::Identity();
	globalPos_ = Eigen::Vector3d::Zero();;

	intrinsicMat_ = Eigen::Matrix3d::Identity();

	distortionParams_.resize(5);
	for (size_t p = 0; p < distortionParams_.size(); ++p)
		distortionParams_[p] = 0;


}

//CameraModel::CameraModel(const std::string& calib_filepath)
//{
//	this->loadCalibrationFile(calib_filepath);
////	this->setCameraType(camera_type);
//}

CameraModel::CameraModel(const CameraModel& copy)
{
	this->intrinsicMat_ = copy.intrinsicMat_;
	this->distortionParams_ = copy.distortionParams_;

	this->intrinsicMat_cv_ = copy.intrinsicMat_cv_.clone();
	this->distortionParams_cv_ = copy.distortionParams_cv_.clone();

	this->globalOrient_mat_ = copy.globalOrient_mat_;
	this->globalPos_ = copy.globalPos_;

	//copy.undistPixToDistPix_x_.copyTo(this->undistPixToDistPix_x_);
	//copy.undistPixToDistPix_y_.copyTo(this->undistPixToDistPix_y_);

	this->cameraName_ = copy.cameraName_;
	this->cameraType_ = copy.cameraType_;
	this->cameraID_ = copy.cameraID_;
	this->imageWidth_ = copy.imageWidth_;
	this->imageHeight_ = copy.imageHeight_;
}

//bool CameraModel::loadCalibrationFile(const std::string& calib_filepath)
//{
//	cv::Mat_<double> temp_orient_vec;
//	cv::Mat_<double> temp_pos_vec;
//	string temp_camera_type_label;
//
//	bool ret = LoadCalibrationFile(calib_filepath, cameraName_, cameraID_, temp_camera_type_label, imageWidth_, imageHeight_, intrinsicMat_, temp_orient_vec, temp_pos_vec, distortionParams_);
//
//	cameraType_ = GetCameraTypeFromLabels(temp_camera_type_label);
//
//	cout << "camera name: " << cameraName_ << endl;
//	cout << "camera ID: " << cameraID_ << endl;
//	cout << "camera type: " << cameraType_ << " - " << temp_camera_type_label << endl;
//
//	cout << "image width: " << imageWidth_ << endl;
//	cout << "image height: " << imageHeight_ << endl;
//
//	// conversion from Nick's representation to general representation
//	cv::Mat_<double> temp_orient_mat;
//	cv::Rodrigues(temp_orient_vec, temp_orient_mat);
//	globalOrient_mat_ = temp_orient_mat.t();
//
//	globalPos_ = -1. * temp_orient_mat * temp_pos_vec;
//
//	cout << "global pos: " << endl << globalPos_ << endl;
////	cout << "global orient: " << endl << globalOrient_vec_ << endl;
//	cout << "global orient: " << endl << globalOrient_mat_ << endl;
//
//	cout << "extrinsics: " << endl << this->extrinsicMat() << endl << endl;
//
//	this->getUndistCoords2DistCoords_();
//
//	return ret;
//}

bool CameraModel::loadCalibrationFile(const std::string& calib_filepath)
{
	ifstream in_calib(calib_filepath.c_str());
	if (in_calib.fail())
	{
		cout << "Error: cannot load calibration file." << endl;
		cout << "\t file path: " << calib_filepath << endl;
		cout << "\t FILE: " << __FILE__ << endl;
		cout << "\t LINE: " << __LINE__ <<  endl;

		return false;
			
	}

	std::string word_buf;
	while (!in_calib.eof())
	{
		word_buf.clear();
		in_calib >> word_buf;

		if (word_buf.empty()) continue;
		else
		{
			if (word_buf == "CameraType")
			{
				//int camera_type;
				//in_calib >> camera_type;
				//this->cameraType_ = (CameraType)camera_type;
				in_calib >> this->cameraType_;
			}
			else if (word_buf == "GlobalPosition")
			{
				in_calib >> this->globalPos_[0] >> this->globalPos_[1] >> this->globalPos_[2];
			}
			else if (word_buf == "GlobalOrient")
			{
				in_calib
					>> this->globalOrient_mat_(0, 0)
					>> this->globalOrient_mat_(0, 1)
					>> this->globalOrient_mat_(0, 2)
					>> this->globalOrient_mat_(1, 0)
					>> this->globalOrient_mat_(1, 1)
					>> this->globalOrient_mat_(1, 2)
					>> this->globalOrient_mat_(2, 0)
					>> this->globalOrient_mat_(2, 1)
					>> this->globalOrient_mat_(2, 2);
			}
			else if (word_buf == "ImageWidth")
			{
				in_calib >> this->imageWidth_ ;
			}
			else if (word_buf == "ImageHeight")
			{
				in_calib >> this->imageHeight_;
			}
			else if (word_buf == "FocalLengthX")
			{
				in_calib >> this->intrinsicMat_(0, 0);
			}
			else if (word_buf == "FocalLengthY")
			{
				in_calib >> this->intrinsicMat_(1, 1);
			}
			else if (word_buf == "PrincipalPointX")
			{
				in_calib >> this->intrinsicMat_(0, 2);
			}
			else if (word_buf == "PrincipalPointY")
			{
				in_calib >> this->intrinsicMat_(1, 2);
			}
			else if (word_buf == "DistortionParams")
			{
				in_calib
					>> this->distortionParams_[0]
					>> this->distortionParams_[1]
					>> this->distortionParams_[2]
					>> this->distortionParams_[3]
					>> this->distortionParams_[4];
			}
			else
			{
				cout << "WTF?" << endl;
			}
		}
	}

	// update opencv version of intrinsic mat and distortion mat
	intrinsicMat_cv_.create(3, 3, CV_64FC1);
	for (size_t r = 0; r < 3; ++r)
	{
		for (size_t c = 0; c < 3; ++c)
		{
			intrinsicMat_cv_.at<double>(r, c) = intrinsicMat_(r, c);
		}
	}

	// distortion params
	distortionParams_cv_.create(5, 1, CV_64FC1);
	for (size_t r = 0; r < 5; ++r)
	{
		distortionParams_cv_.at<double>(r, 0) = distortionParams_(r, 0);
	}

	return true;
}

bool CameraModel::writeCalibrationFile(const std::string& calib_filepath) const
{
	ofstream out_txt(calib_filepath.c_str());
	if (out_txt.fail())
	{
		cout << "Error: cannot write calibraion file." << endl;
		cout << "\t file path: " << calib_filepath << endl;
		cout << "\t FILE: " << __FILE__ << endl;
		cout << "\t LINE: " << __LINE__ << endl;
		return false;
	}

	// meta information
	out_txt << "CameraType" << " " << this->cameraType_ << endl;

	// extrinsics
	out_txt << "GlobalPosition" << " " << this->globalPos_[0] << " " << this->globalPos_[1] << " " << this->globalPos_[2] << endl;
	out_txt << "GlobalOrient" << " "
		<< this->globalOrient_mat_(0, 0) << " "
		<< this->globalOrient_mat_(0, 1) << " "
		<< this->globalOrient_mat_(0, 2) << " "
		<< this->globalOrient_mat_(1, 0) << " "
		<< this->globalOrient_mat_(1, 1) << " "
		<< this->globalOrient_mat_(1, 2) << " "
		<< this->globalOrient_mat_(2, 0) << " "
		<< this->globalOrient_mat_(2, 1) << " "
		<< this->globalOrient_mat_(2, 2) << endl;

	// intrinsics
	out_txt << "ImageWidth" << " " << this->imageWidth_ << endl;
	out_txt << "ImageHeight" << " " << this->imageHeight_ << endl;
	out_txt << "FocalLengthX" << " " << this->intrinsicMat_(0, 0) << endl;
	out_txt << "FocalLengthY" << " " << this->intrinsicMat_(1, 1) << endl;
	out_txt << "PrincipalPointX" << " " << this->intrinsicMat_(0, 2) << endl;
	out_txt << "PrincipalPointY" << " " << this->intrinsicMat_(1, 2) << endl;
	out_txt << "DistortionParams" << " "
		<< this->distortionParams_[0] << " "
		<< this->distortionParams_[1] << " "
		<< this->distortionParams_[2] << " "
		<< this->distortionParams_[3] << " "
		<< this->distortionParams_[4] << endl;


	return true;
}

void CameraModel::setImageResolution(size_t image_width, size_t image_height)
{
	imageWidth_ = image_width;
	imageHeight_ = image_height;
}

void CameraModel::setIntrinsicParams(double focal_x, double focal_y, double pp_x, double pp_y)
{
	intrinsicMat_ = Eigen::Matrix3d::Identity();
	intrinsicMat_(0, 0) = focal_x;
	intrinsicMat_(1, 1) = focal_y;
	intrinsicMat_(0, 2) = pp_x;
	intrinsicMat_(1, 2) = pp_y;

	intrinsicMat_cv_.create(3, 3, CV_64FC1);
	for (size_t r = 0; r < 3; ++r)
	{
		for (size_t c = 0; c < 3; ++c)
		{
			intrinsicMat_cv_.at<double>(r, c) = intrinsicMat_(r, c);
		}
	}
}

void CameraModel::setIntrinsicMat(const Eigen::Matrix3d& intrinsic_mat)
{
	intrinsicMat_ = intrinsic_mat;

	intrinsicMat_cv_.create(3, 3, CV_64FC1);
	for (size_t r = 0; r < 3; ++r)
	{
		for (size_t c = 0; c < 3; ++c)
		{
			intrinsicMat_cv_.at<double>(r, c) = intrinsicMat_(r, c);
		}
	}

}

Eigen::MatrixXd CameraModel::extrinsicMat() const
{
	Eigen::MatrixXd extrinsic_mat(3, 4);
	Eigen::Vector3d local_t = -1. * globalOrient_mat_ * globalPos_;

	for (size_t r = 0; r < 3; ++r)
	{
		for (size_t c = 0; c < 3; ++c)
		{
			extrinsic_mat(r, c) = globalOrient_mat_(r, c);
		}
		extrinsic_mat(r, 3) = local_t(r, 0);
	}
	return extrinsic_mat;
}

double CameraModel::computeDepth(const Eigen::Vector3d& global_pos_3D) const
{
	return (this->convertGlobalPosToCameraLocalPos(global_pos_3D))[0];
}

cv::Mat CameraModel::extrinsicMat_cv() const
{
	Eigen::MatrixXd ext_mat_eigen = this->extrinsicMat();

	cv::Mat ext_mat_cv(3, 4, CV_64FC1);
	for (size_t r = 0; r < 3; ++r)
	{
		for (size_t c = 0; c < 4; ++c)
		{
			ext_mat_cv.at<double>(r, c) = ext_mat_eigen(r, c);
		}
	}

	return ext_mat_cv;

}


void CameraModel::setExtrinsicMat(const Eigen::MatrixXd& extrinsic_mat)
{
	for (size_t r = 0; r < 3; ++r)
		for (size_t c = 0; c < 3; ++c)
			globalOrient_mat_(r, c) = extrinsic_mat(r, c);

	Eigen::Vector3d local_t = extrinsic_mat.col(3);
	
	globalPos_ = -1. * globalOrient_mat_.transpose() * local_t;
}

void CameraModel::setExtrinsicParams(const Eigen::Matrix3d& global_orient, const Eigen::Vector3d& global_pos)
{
	globalOrient_mat_ = global_orient;
	globalPos_ = global_pos;
}

/*
void CameraModel::getUndistCoords2DistCoords_()
{
	cv::initUndistortRectifyMap(intrinsicMat_, distortionParams_, cv::Mat(), intrinsicMat_, cv::Size(imageWidth_, imageHeight_), CV_32FC1, undistPixToDistPix_x_, undistPixToDistPix_y_);


	//// debug
	//double k1 = distortionParams_(0, 0);
	//double k2 = distortionParams_(1, 0);
	//double p1 = distortionParams_(2, 0);
	//double p2 = distortionParams_(3, 0);

	//int target_undist_x = 102;
	//int target_undist_y = 134;
	//cv::Mat_<double> undist_pix_pos(3, 1);
	//undist_pix_pos(0, 0) = target_undist_x;
	//undist_pix_pos(1, 0) = target_undist_y;
	//undist_pix_pos(2, 0) = 1;

	//cv::Mat_<double> undist_norm_pos(3, 1);
	//undist_norm_pos = intrinsicMat_.inv() * undist_pix_pos;

	//double r2 = undist_norm_pos(0, 0) * undist_norm_pos(0, 0) + undist_norm_pos(1, 0) * undist_norm_pos(1, 0);

	//double curr_x = undist_norm_pos(0, 0);
	//double curr_y = undist_norm_pos(1, 0);
	//cv::Mat_<double> dist_norm_pos(3, 1);
	//dist_norm_pos(0, 0) = curr_x * (1 + k1 * r2 + k2 * r2 * r2)
	//	+ 2. * p1 * curr_x * curr_y + p2 * (r2 + 2 * curr_x * curr_x);
	//dist_norm_pos(1, 0) = curr_y * (1 + k1 * r2 + k2 * r2 * r2)
	//	+ p1 * (r2 + 2 * curr_y * curr_y) + 2. * p2 * curr_x * curr_y;
	//dist_norm_pos(2, 0) = 1;

	//cv::Mat_<double> dist_pix_pos(3, 1);
	//dist_pix_pos = intrinsicMat_ * dist_norm_pos;

	//cout << "dist x: " << dist_pix_pos(0, 0) << endl;
	//cout << "dist y: " << dist_pix_pos(1, 0) << endl;

	//cout << "in map..." << endl;
	//cout << "dist x: " << undistPixToDistPix_x_(target_undist_y, target_undist_x) << endl;
	//cout << "dist y: " << undistPixToDistPix_y_(target_undist_y, target_undist_x) << endl;
}

cv::Point2d CameraModel::getDistPixPosFromUndistPixPos(const cv::Point2d& undistort_pos_pix) const
{
	cv::Point2d ret_dist_pos;
	ret_dist_pos.x = undistPixToDistPix_x_(undistort_pos_pix.y, undistort_pos_pix.x);
	ret_dist_pos.y = undistPixToDistPix_y_(undistort_pos_pix.y, undistort_pos_pix.x);

	return ret_dist_pos;
}
*/

//
//bool CameraModel::applyUndistortion()
//{
//	//if (!(this->hasImage()))
//	//	return false;
//
//	cv::Mat undistort_image(image_.rows, image_.cols, image_.type());
//	cv::undistort( image_, undistort_image, intrinsicMat_, distortionParams_);
//
//	undistort_image.copyTo(undistort_image);
//
//	return true;
//}

//////////////////////////////////////
// for color camera
ColorCameraModel::ColorCameraModel()
	: CameraModel()
{
	this->cameraType_ = Camera_RGB;
}

ColorCameraModel::ColorCameraModel(const ColorCameraModel& copy)
	: CameraModel(copy)
{
	this->cameraType_ = Camera_RGB;
	copy.colorImage_.copyTo(this->colorImage_);
}

ColorCameraModel::~ColorCameraModel()
{
}

bool ColorCameraModel::loadImage(const std::string& image_filepath)
{
	colorImage_ = cv::imread(image_filepath);

	return true;
}

bool ColorCameraModel::writeImage(const std::string& image_filepath) const
{
	cv::imwrite(image_filepath, colorImage_);

	return true;
}

//bool ColorCameraModel::loadImageSeqList(const std::string& directory, const std::string& image_list_filename)
//{
//	vector<string> image_filenames;
//	bool file_loaded = LoadListFile(directory + "/" + image_list_filename, image_filenames);
//	if (!file_loaded)
//	{
//		cout << "Error in loading image list." << endl;
//		cout << "\t directory: " << directory << endl;
//		cout << "\t list file name: " << image_list_filename << endl;
//		cout << "\t FILE: " << __FILE__ << endl;
//		cout << "\t LINE: " << __LINE__ << endl;
//
//		return false;
//	}
//
//
//	imageFilePathSeq_.clear();
//	for (size_t i = 0; i < image_filenames.size(); ++i)
//	{
//		imageFilePathSeq_.push_back( directory + "/" + image_filenames[i] );
//	}
//
//
//	return true;
//}

//////////////////////////////////////
// for depth camera

DepthCameraModel::DepthCameraModel()
: CameraModel(), minDepth_m_(0.01), maxDepth_m_( 15 )
{
	this->cameraType_ = Camera_IRDEPTH;
}

/*
DepthCameraModel::DepthCameraModel(const string& calibration_filepath)
: CameraModel(calibration_filepath), minDepth_m_(0.01), maxDepth_m_(15)
{
}
*/

DepthCameraModel::DepthCameraModel(const DepthCameraModel& copy)
: CameraModel(copy), minDepth_m_(0.01), maxDepth_m_(15)
{
	copy.depthImage_m_.copyTo(this->depthImage_m_);
//	copy.normalMap_.copyTo(this->normalMap_);
}

bool DepthCameraModel::loadImage(const string& png_filepath)
{
	//LoadDepthMapPNG(png_filepath, depthImage_m_, MiliMeter2Meter);
	cout << "Warning: DepthCameraModel::loadImage() is not implemented yet..." << endl;
	cout << "\t FILE: " << __FILE__ << endl;
	cout << "\t LINE: " << __LINE__ << endl;
	return false;
}

//void DepthCameraModel::createValidDepthMask_()
//{
//	validDepthMask_.create(depthImage_m_.rows, depthImage_m_.cols);
//
//	cout << "Now checking depth validity...";
//	for (size_t r = 0; r < depthImage_m_.rows; ++r)
//	{
//		for (size_t c = 0; c < depthImage_m_.cols; ++c)
//		{
//			if (minDepth_m_ < depthImage_m_(r, c) && depthImage_m_(r, c) < maxDepth_m_)
//				validDepthMask_(r, c) = true;
//			else
//				validDepthMask_(r, c) = false;
//		}
//	}
//	cout << "done!" << endl;
//
//}



bool DepthCameraModel::writePseudoDepthImage(const string& png_filepath, const double scale ) const
{

//	cv::Mat_<float> temp_depth_image = scale * depthImage_m_;
	cv::Mat pseudo_depth_image(depthImage_m_.rows, depthImage_m_.cols, CV_8UC1);
	for (size_t r = 0; r < pseudo_depth_image.rows; ++r)
	{
		for (size_t c = 0; c < pseudo_depth_image.cols; ++c)
		{
			pseudo_depth_image.at<uchar>(r, c) = (uchar)(scale * depthImage_m_(r,c));
		}
	}

	cv::imwrite(png_filepath, pseudo_depth_image);

	return true;
}

/*
void DepthCameraModel::applySmoothing()
{
	cv::Mat_<float> filtered_image(depthImage_m_.rows, depthImage_m_.cols);

	cv::bilateralFilter(depthImage_m_, filtered_image, 2, .1, 3);

	filtered_image.copyTo(depthImage_m_);
}

static int debug_id = 0;

void DepthCameraModel::computeNormalMap()
{
//	normalMap_.create(depthImage_m_.rows, depthImage_m_.cols);
	//normalMap_.create(3, imageWidth_ * imageHeight_);
	normalMap_.create(imageWidth_ * imageHeight_, 3);

	cv::Mat_<double> point_pos_3D;
	this->convertDepthMapToPointCloud(point_pos_3D);


	for (size_t curr_y = 0; curr_y < imageHeight_; ++curr_y)
	{
		for (size_t curr_x = 0; curr_x < imageWidth_; ++curr_x)
		{
			if (minDepth_m_ < depthImage_m_(curr_y, curr_x) && depthImage_m_(curr_y, curr_x) < maxDepth_m_)
			{

				int next_x = curr_x + 1;
				if (curr_x == imageWidth_ - 1) next_x = imageWidth_ - 2;

				int next_y = curr_y + 1;
				if (curr_y == imageHeight_ - 1) next_y = imageHeight_ - 2;

				double curr_depth = depthImage_m_(curr_y, curr_x);
				double next_depth1 = depthImage_m_(curr_y, next_x);
				double next_depth2 = depthImage_m_(next_y, curr_x);

				static const double depth_diff_thresh = 0.2;
				if (hasDepthInRange_(next_depth1) && hasDepthInRange_(next_depth2)
					&& fabs(curr_depth - next_depth1) < depth_diff_thresh
					&& fabs(curr_depth - next_depth2) < depth_diff_thresh
					&& fabs(next_depth2 - next_depth1) < depth_diff_thresh)
				{
					size_t target_pt_idx;
					target_pt_idx = array1DIdxFromPixXY_(curr_x, curr_y); //curr_y * depthImage_m_.cols + curr_x;
					cv::Mat_<double> curr_pos(3, 1);
					curr_pos(0, 0) = point_pos_3D(target_pt_idx, 0);
					curr_pos(1, 0) = point_pos_3D(target_pt_idx, 1);
					curr_pos(2, 0) = point_pos_3D(target_pt_idx, 2);
					//curr_pos(0, 0) = point_pos_3D(0, target_pt_idx);
					//curr_pos(1, 0) = point_pos_3D(1, target_pt_idx);
					//curr_pos(2, 0) = point_pos_3D(2, target_pt_idx);

					target_pt_idx = array1DIdxFromPixXY_(next_x, curr_y);//curr_y * depthImage_m_.cols + next_x;
					cv::Mat_<double> next_pos1(3, 1);
					next_pos1(0, 0) = point_pos_3D(target_pt_idx, 0);
					next_pos1(1, 0) = point_pos_3D(target_pt_idx, 1);
					next_pos1(2, 0) = point_pos_3D(target_pt_idx, 2);
					//next_pos1(0, 0) = point_pos_3D(0, target_pt_idx);
					//next_pos1(1, 0) = point_pos_3D(1, target_pt_idx);
					//next_pos1(2, 0) = point_pos_3D(2, target_pt_idx);

					target_pt_idx = array1DIdxFromPixXY_(curr_x, next_y);//next_y * depthImage_m_.cols + curr_x;
					cv::Mat_<double> next_pos2(3, 1);
					next_pos2(0, 0) = point_pos_3D(target_pt_idx, 0);
					next_pos2(1, 0) = point_pos_3D(target_pt_idx, 1);
					next_pos2(2, 0) = point_pos_3D(target_pt_idx, 2);
					//next_pos2(0, 0) = point_pos_3D(0, target_pt_idx);
					//next_pos2(1, 0) = point_pos_3D(1, target_pt_idx);
					//next_pos2(2, 0) = point_pos_3D(2, target_pt_idx);

					cv::Mat_<double> vec_side = next_pos1 - curr_pos;
					cv::Mat_<double> vec_vert = next_pos2 - curr_pos;
					cv::Mat_<double> normal = vec_vert.cross(vec_side);
					normal /= (double)(sqrt(normal.dot(normal)));


					cv::Mat_<double> pt_to_cam = curr_pos - this->globalPos_;

					if (normal.dot(pt_to_cam) < 0)
						normal *= -1.;

//					normalMap_(curr_y, curr_x) = cv::Vec3d(normal(0, 0), normal(1, 0), normal(2, 0));
					size_t row_idx = array1DIdxFromPixXY_(curr_x, curr_y);
					normalMap_(row_idx,0) = normal(0, 0);
					normalMap_(row_idx,1) = normal(1, 0);
					normalMap_(row_idx,2) = normal(2, 0);
					//normalMap_(0, col_idx) = normal(0, 0);
					//normalMap_(1, col_idx) = normal(1, 0);
					//normalMap_(2, col_idx) = normal(2, 0);

				}
			}
			else
			{
//				normalMap_(curr_y, curr_x) = cv::Vec3d(0, 0, 0);
				size_t row_idx = array1DIdxFromPixXY_(curr_x, curr_y);
				normalMap_(row_idx,0) = 0;
				normalMap_(row_idx,1) = 0;
				normalMap_(row_idx,2) = 0;
			}


		}
	}

////	cv::Mat temp_normal_map(depthImage_m_.rows, depthImage_m_.cols, CV_8UC3);
//	for (size_t y = 0; y < imageHeight_; ++y)
//	{
//		for (size_t x = 0; x < imageWidth_; ++x)
//		{
//			if (250 < x && x < 300
//				&& 160 < y && y < 180)
//			{
//				double temp_R = (normalMap_(y, x)[0]);
//				double temp_G = (normalMap_(y, x)[1]);
//				double temp_B = (normalMap_(y, x)[2]);
//				cout << "r,g,b = " << temp_R << " " << temp_G << " " << temp_B << endl;
//
//				//cv::Mat_<double> global_normal(3, 1);
//				//global_normal(0, 0) = normalMap_(y, x)[0];
//				//global_normal(1, 0) = normalMap_(y, x)[1];
//				//global_normal(2, 0) = normalMap_(y, x)[2];
//				//cv::Mat_<double> local_normal = globalOrient_mat_ * global_normal;
//
//				//double normal_x = local_normal(0, 0);
//				//double normal_y = local_normal(1, 0);
//				//double normal_z = local_normal(2, 0);
//
//				//if ((normal_x * normal_x + normal_y * normal_y + normal_z* normal_z) > 0.001)
//				//{
//
//				//	temp_normal_map.at<cv::Vec3b>(y, x)[0] = (int)(255 - normal_z * 255);
//				//	temp_normal_map.at<cv::Vec3b>(y, x)[1] = (int)(255 - normal_y * 255);
//				//	temp_normal_map.at<cv::Vec3b>(y, x)[2] = (int)(255 - normal_x * 255);
//				//}
//			}
//		}
//	}
//
//	normalMap_ *= 255.;
//	string debug_normal_fpath = "normal_" + to_string(debug_id) + ".png";
//	cout << "debug_normal_fpath: " << debug_normal_fpath << endl;
//	cv::imwrite(debug_normal_fpath, normalMap_);
//	++debug_id;
}


bool DepthCameraModel::applyUndistortion()
{

	cv::Mat_<float> undist_depth_map(depthImage_m_.rows, depthImage_m_.cols);
	for (size_t undist_y = 0; undist_y < undist_depth_map.rows; ++undist_y)
	{
		for (size_t undist_x = 0; undist_x < depthImage_m_.cols; ++undist_x)
		{
			//int dist_x = (int)(round(undistPixToDistPix_x_(undist_y, undist_x)));
			//int dist_y = (int)(round(undistPixToDistPix_y_(undist_y, undist_x)));
			int dist_x = (int)((undistPixToDistPix_x_(undist_y, undist_x)));
			int dist_y = (int)((undistPixToDistPix_y_(undist_y, undist_x)));

			// nearest neighbor, not bilinear interpoation
			undist_depth_map(undist_y, undist_x) = depthImage_m_(dist_y, dist_x);
		}
	}

	undist_depth_map.copyTo(depthImage_m_);

	return true;
}

void DepthCameraModel::convertDepthMapToPointCloud(cv::Mat_<double>& point_pos)//std::vector<cv::Point3f>& point_pos)
{
	point_pos.create(depthImage_m_.rows * depthImage_m_.cols, 3);
	//point_pos.create(3, depthImage_m_.rows * depthImage_m_.cols);

	cv::Mat_<double> inv_intrinsics = intrinsicMat_.inv();

	for (size_t undist_y = 0; undist_y < imageHeight_; ++undist_y)
	{
		for (size_t undist_x = 0; undist_x < imageWidth_; ++undist_x)
		{
			double curr_depth = depthImage_m_(undist_y, undist_x);


			if (curr_depth <= minDepth_m_ || curr_depth >= maxDepth_m_ )
				curr_depth = 0;


			{
				cv::Mat_<double> pixel_pos(3, 1);
				pixel_pos(0, 0) = undist_x;
				pixel_pos(1, 0) = undist_y;
				pixel_pos(2, 0) = 1;

				cv::Mat_<double> normalized_pixel_pos(3, 1);
				normalized_pixel_pos = inv_intrinsics * pixel_pos;

				cv::Mat_<double> point_pos_in_cam(3, 1);
				point_pos_in_cam = curr_depth * normalized_pixel_pos;


				cv::Mat_<double> global_point_pos;

				global_point_pos = globalOrient_mat_.t() * point_pos_in_cam + globalPos_;

				size_t curr_point_idx = undist_y * depthImage_m_.cols + undist_x;
				//size_t curr_point_idx = this->colIdxFromPixXY_(undist_x, undist_y);
				point_pos(curr_point_idx, 0) = global_point_pos(0, 0);
				point_pos(curr_point_idx, 1) = global_point_pos(1, 0);
				point_pos(curr_point_idx, 2) = global_point_pos(2, 0);
				//point_pos(0, curr_point_idx) = global_point_pos(0, 0);
				//point_pos(1, curr_point_idx) = global_point_pos(1, 0);
				//point_pos(2, curr_point_idx) = global_point_pos(2, 0);

				//cout << curr_point_idx << endl;
				//cout << "global point pos " << global_point_pos << endl;
			}
		}
	}

}

cv::Mat_<double> DepthCameraModel::computePointPosAt(size_t x, size_t y) const
{

	cv::Mat_<double> normalized_pos;
	normalized_pos(0,0) = (x - intrinsicMat_(0, 2)) / intrinsicMat_(0, 0);
	normalized_pos(1,0) = (y - intrinsicMat_(1, 2)) / intrinsicMat_(1, 1);
	normalized_pos(2,0) = 1;

	const double curr_depth = depthImage_m_(y, x);
	cv::Mat_<double> point_pos_in_cam(3, 1);
	point_pos_in_cam = curr_depth * normalized_pos;

	cv::Mat_<double> global_point_pos(3, 1);
	global_point_pos = globalOrient_mat_.t() * point_pos_in_cam + globalPos_;

	return global_point_pos;
 }
 */

//bool DepthCameraModel::loadDepthImageSeqList(const std::string& directory, const std::string& image_list_filename)
//{
//	vector<string> image_filenames;
//	bool file_loaded = LoadListFile(directory + "/" + image_list_filename, image_filenames);
//	if (!file_loaded)
//	{
//		cout << "Error in loading depth image list." << endl;
//		cout << "\t directory: " << directory << endl;
//		cout << "\t list file name: " << image_list_filename << endl;
//		cout << "\t FILE: " << __FILE__ << endl;
//		cout << "\t LINE: " << __LINE__ << endl;
//
//		return false;
//	}
//
//
//	depthImageFilePathSeq_.clear();
//	for (size_t i = 0; i < image_filenames.size(); ++i)
//	{
//		depthImageFilePathSeq_.push_back(directory + "/" + image_filenames[i]);
//	}
//
//
//	return true;
//}
//
//bool DepthCameraModel::loadIRImageSeqList(const std::string& directory, const std::string& image_list_filename)
//{
//	vector<string> image_filenames;
//	bool file_loaded = LoadListFile(directory + "/" + image_list_filename, image_filenames);
//	if (!file_loaded)
//	{
//		cout << "Error in loading depth image list." << endl;
//		cout << "\t directory: " << directory << endl;
//		cout << "\t list file name: " << image_list_filename << endl;
//		cout << "\t FILE: " << __FILE__ << endl;
//		cout << "\t LINE: " << __LINE__ << endl;
//
//		return false;
//	}
//
//
//	irImageFilePathSeq_.clear();
//	for (size_t i = 0; i < image_filenames.size(); ++i)
//	{
//		irImageFilePathSeq_.push_back(directory + "/" + image_filenames[i]);
//	}
//
//
//	return true;
//}


Eigen::Vector3d CameraModel::globalRayDirOfPixel(const Eigen::Vector2d& pixel_pos) const
{
	//cout << "WARNING: this needs to consider distortion params!" << endl;
	//Eigen::Vector3d local_ray = intrinsicMat_.inverse() * pixel_pos.homogeneous();
	//Eigen::Vector3d global_ray = globalOrient_mat_.transpose() * local_ray;
	//global_ray.normalize();
	//return global_ray;

	// pixel location
	//cv::Mat_<double> pixel_pos_cv(2, 1);
	//pixel_pos_cv(0, 0) = pixel_pos[0];
	//pixel_pos_cv(1, 0) = pixel_pos[1];
	cv::Mat pixel_pos_cv(1, 1, CV_64FC2);
	pixel_pos_cv.at<cv::Vec2d>(0, 0)[0] = pixel_pos[0];
	pixel_pos_cv.at<cv::Vec2d>(0, 0)[1] = pixel_pos[1];

	// undistortion and normalization 
	//cv::Mat_<double> undist_norm_pos_cv(2, 1);

	//cout << "intrinsic rows: " << intrinsicMat_cv_.rows << endl;
	//cout << "intrinsic cols: " << intrinsicMat_cv_.cols << endl;
	//cout << "intrinsic mat: " << endl << intrinsicMat_cv_ << endl;

	cv::Mat undist_norm_pos_cv(1, 1, CV_64FC2);
	cv::undistortPoints(pixel_pos_cv, undist_norm_pos_cv, intrinsicMat_cv_, distortionParams_cv_, cv::noArray(), cv::noArray());

	// back to Eigen
	Eigen::Vector3d local_ray;
	//local_ray[0] = undist_norm_pos_cv(0, 0);
	//local_ray[1] = undist_norm_pos_cv(1, 0);
	local_ray[0] = undist_norm_pos_cv.at<cv::Vec2d>(0, 0)[0];
	local_ray[1] = undist_norm_pos_cv.at<cv::Vec2d>(0, 0)[1];
	local_ray[2] = 1;

	Eigen::Vector3d global_ray = globalOrient_mat_.transpose() * local_ray;
	global_ray.normalize();
	return global_ray;
}

Eigen::Vector3d CameraModel::undistortPixelPos_norm(const Eigen::Vector3d& dist_pos_norm) const
{
	cv::Mat dist_pos_norm_cv(1, 1, CV_64FC2);
	dist_pos_norm_cv.at<cv::Vec2d>(0, 0)[0] = dist_pos_norm[0];
	dist_pos_norm_cv.at<cv::Vec2d>(0, 0)[1] = dist_pos_norm[1];

	cv::Mat undist_pos_norm_cv(1, 1, CV_64FC2);
	cv::undistortPoints(dist_pos_norm_cv, undist_pos_norm_cv, cv::Mat::eye(3,3,CV_64FC1), distortionParams_cv_, cv::noArray(), cv::noArray());

	return Eigen::Vector3d(undist_pos_norm_cv.at<cv::Vec2d>(0, 0)[0], undist_pos_norm_cv.at<cv::Vec2d>(0, 0)[1], 1);
}

Eigen::Vector3d CameraModel::undistortPixelPos_pix(const Eigen::Vector3d& dist_pos_pix) const
{
	Eigen::Vector3d dist_pos_norm = intrinsicMat_.inverse() * dist_pos_pix;

	cv::Mat dist_pos_norm_cv(1, 1, CV_64FC2);
	dist_pos_norm_cv.at<cv::Vec2d>(0, 0)[0] = dist_pos_norm[0];
	dist_pos_norm_cv.at<cv::Vec2d>(0, 0)[1] = dist_pos_norm[1];

	cv::Mat undist_pos_norm_cv(1, 1, CV_64FC2);
	cv::undistortPoints(dist_pos_norm_cv, undist_pos_norm_cv, cv::Mat::eye(3, 3, CV_64FC1), distortionParams_cv_, cv::noArray(), cv::noArray());

	Eigen::Vector3d undist_pos_norm
		= Eigen::Vector3d(undist_pos_norm_cv.at<cv::Vec2d>(0, 0)[0], undist_pos_norm_cv.at<cv::Vec2d>(0, 0)[1], 1);

	return intrinsicMat_ * undist_pos_norm;
}


Eigen::Vector3d CameraModel::globalRayDirOfPixel_NoDist(const Eigen::Vector2d& pixel_pos) const
{
	//cout << "WARNING: this needs to consider distortion params!" << endl;
	Eigen::Vector3d local_ray = intrinsicMat_.inverse() * pixel_pos.homogeneous();
	Eigen::Vector3d global_ray = globalOrient_mat_.transpose() * local_ray;
	global_ray.normalize();
	return global_ray;

}


Eigen::Vector3d CameraModel::backProjectTo3DWithDistance(const Eigen::Vector2d& pixel_pos, const double distance) const
{
	cout << "WARNING: this needs to consider distortion params!" << endl;

	Eigen::Vector3d local_ray = intrinsicMat_.inverse() * pixel_pos.homogeneous();
	Eigen::Vector3d global_ray = globalOrient_mat_.transpose() * local_ray;
	global_ray.normalize();

	return globalPos_ + distance * global_ray;
}

Eigen::Vector3d CameraModel::computeGlobalPosFromDepth_NoDist(const Eigen::Vector2d& pixel_pos, const double depth) const
{
	Eigen::Vector3d local_pos = depth * intrinsicMat_.inverse() * pixel_pos.homogeneous();
	
	return globalPos_ + globalOrient_mat_.transpose() * local_pos;
}

Eigen::Vector3d CameraModel::computeGlobalPosFromDepth_NoDist(const Eigen::Vector3d& pixel_pos_homo, const double depth) const
{
	Eigen::Vector3d local_pos = depth * intrinsicMat_.inverse() * pixel_pos_homo;

	return globalPos_ + globalOrient_mat_.transpose() * local_pos;
}


#include "PerspectiveProjection.h"

Eigen::Vector3d CameraModel::convertGlobalPosToCameraLocalPos(const Eigen::Vector3d& global_pos_3D) const
{
	return TransformGlobalPosIntoCameraCoordSystem(global_pos_3D, this->extrinsicMat());
}

Eigen::Vector3d CameraModel::projectPoint3DToImageNoDist_pix(const Eigen::Vector3d& pos_3D) const
{
	//Eigen::Vector3d pos_pix_homo = TransformGlobalPos3DToPixelCoords(
	//	pos_3D,
	//	this->extrinsicMat(),
	//	this->intrinsicMat_);
//	return Eigen::Vector2d(pos_pix_homo[0], pos_pix_homo[1]);

	return TransformGlobalPos3DToPixelCoords(
		pos_3D,
		this->extrinsicMat(),
		this->intrinsicMat_);

}

Eigen::Vector3d CameraModel::projectPoint3DToImageWithDist_pix(const Eigen::Vector3d& pos_3D) const
{
	Eigen::Vector3d normalized_pos = TransformGlobalPos3DToNormalizedImageCoords(pos_3D, this->extrinsicMat());

	//// apply distortion
	//const double r2 = normalized_pos[0] * normalized_pos[0] + normalized_pos[1] * normalized_pos[1];
	//const double r4 = r2 * r2;
	//const double r6 = r4 * r2;
	//const double k1 = distortionParams_[0];
	//const double k2 = distortionParams_[1];
	//const double p1 = distortionParams_[2];
	//const double p2 = distortionParams_[3];
	//const double k3 = distortionParams_[4];

	//Eigen::Vector2d dist_pos_pix;
	//dist_pos_pix[0]
	//	= normalized_pos[0] * (1 + k1 * r2 + k2 * r4 + k3 * r6)
	//	+ 2. * p1 * normalized_pos[0] * normalized_pos[1] + p2 * (r2 + 2. * normalized_pos[0] * normalized_pos[0]);

	//dist_pos_pix[1]
	//	= normalized_pos[1] * (1 + k1 * r2 + k2 * r4 + k3 * r6)
	//	+ p1 * (r2 + 2. * normalized_pos[1] * normalized_pos[1]) + 2. * p2 * normalized_pos[0] * normalized_pos[1];


	//return dist_pos_pix;

	return ( intrinsicMat_ * this->applyRadialDistortion_norm(normalized_pos));
}

Eigen::Vector3d CameraModel::applyRadialDistortion_norm(const Eigen::Vector3d& undist_pos_norm) const
{
	// apply distortion
	const double r2 = undist_pos_norm[0] * undist_pos_norm[0] + undist_pos_norm[1] * undist_pos_norm[1];
	const double r4 = r2 * r2;
	const double r6 = r4 * r2;
	const double k1 = distortionParams_[0];
	const double k2 = distortionParams_[1];
	const double p1 = distortionParams_[2];
	const double p2 = distortionParams_[3];
	const double k3 = distortionParams_[4];


	Eigen::Vector3d dist_pos_norm;
	dist_pos_norm[0]
		= undist_pos_norm[0] * (1 + k1 * r2 + k2 * r4 + k3 * r6)
		+ 2. * p1 * undist_pos_norm[0] * undist_pos_norm[1] + p2 * (r2 + 2. * undist_pos_norm[0] * undist_pos_norm[0]);

	dist_pos_norm[1]
		= undist_pos_norm[1] * (1 + k1 * r2 + k2 * r4 + k3 * r6)
		+ p1 * (r2 + 2. * undist_pos_norm[1] * undist_pos_norm[1]) + 2. * p2 * undist_pos_norm[0] * undist_pos_norm[1];

	dist_pos_norm[2] = 1;

	return dist_pos_norm;

}

bool CameraModel::isPointInFrontOfCamera(const Eigen::Vector3d& pos_3D, const double thresh_local_z) const
{
	Eigen::Vector3d pos_in_cam_coord = TransformGlobalPosIntoCameraCoordSystem(pos_3D, this->extrinsicMat());
	return (pos_in_cam_coord[2] > thresh_local_z) ? true : false;
}

void DepthCameraModel::reconstructPointCloud(const cv::Mat_<float>& depth_map_m, const double max_depth_m, vector<Eigen::Vector3d>& out_point_cloud)
{
	out_point_cloud.clear();
	const Eigen::Matrix3d inv_intrinsics = intrinsicMat_.inverse();

	for (size_t y = 0; y < depth_map_m.rows; ++y)
	{
		for (size_t x = 0; x < depth_map_m.cols; ++x)
		{
			const double curr_depth = (double)(depth_map_m(y, x));

			if (0.001 < curr_depth && curr_depth < max_depth_m)
			{
				Eigen::Vector3d dist_pos_pix(x, y, 1);
				Eigen::Vector3d dist_pos_norm = inv_intrinsics * dist_pos_pix;
				Eigen::Vector3d undist_pos_norm = this->undistortPixelPos_norm(dist_pos_norm);

				Eigen::Vector3d local_pos_3D = undist_pos_norm * curr_depth;

				Eigen::Vector3d global_pos_3D = this->globalOrient_mat_.transpose() * local_pos_3D + this->globalPos_;

				out_point_cloud.push_back(global_pos_3D);

			}
		}
	}

}

void DepthCameraModel::reconstructPointCloudNormals(
	const cv::Mat_<float>& depth_map_m, const double max_depth_m,
	const int boundary_margin_pix,
	vector<Eigen::Vector3d>& out_point_cloud,
	vector<Eigen::Vector3d>& out_normals)
{
	out_point_cloud.clear();
	out_normals.clear();

	const Eigen::Matrix3d inv_intrinsics = intrinsicMat_.inverse();

	// point map initialization
	vector<vector<Eigen::Vector3d>> point_map;
	point_map.resize(depth_map_m.rows);
	for (size_t y = 0; y < depth_map_m.rows; ++y)
	{
		point_map[y].resize(depth_map_m.cols, Eigen::Vector3d::Zero());
	}

	// compute point map
	for (size_t y = boundary_margin_pix; y < depth_map_m.rows - boundary_margin_pix; ++y)
	{
		for (size_t x = boundary_margin_pix; x < depth_map_m.cols - boundary_margin_pix; ++x)
		{
			const double curr_depth = (double)(depth_map_m(y, x));

			if (0.001 < curr_depth && curr_depth < max_depth_m)
			{
				Eigen::Vector3d dist_pos_pix(x, y, 1);
				Eigen::Vector3d dist_pos_norm = inv_intrinsics * dist_pos_pix;
				Eigen::Vector3d undist_pos_norm = this->undistortPixelPos_norm(dist_pos_norm);
				Eigen::Vector3d local_pos_3D = undist_pos_norm * curr_depth;
				Eigen::Vector3d global_pos_3D = this->globalOrient_mat_.transpose() * local_pos_3D + this->globalPos_;

				//				out_point_cloud.push_back(global_pos_3D);

				point_map[y][x] = global_pos_3D;
			}
		}
	}

	// prepare for outputing point cloud and corresponding normals
	for (size_t y = boundary_margin_pix; y < point_map.size() - boundary_margin_pix; ++y)
	{
		for (size_t x = boundary_margin_pix; x < point_map[y].size() - boundary_margin_pix; ++x)
		{
			// if the element is not zero vector...
			if ((point_map[y][x] - Eigen::Vector3d::Zero()).norm() > 0.0001)
			{
				if ((point_map[y + 1][x] - Eigen::Vector3d::Zero()).norm() > 0.0001
					&& (point_map[y][x + 1] - Eigen::Vector3d::Zero()).norm() > 0.0001)
				{
					Eigen::Vector3d dp_dy = point_map[y + 1][x] - point_map[y][x];
					Eigen::Vector3d dp_dx = point_map[y][x + 1] - point_map[y][x];

					if (dp_dy.norm() < 0.2 && dp_dx.norm() < 0.2)
					{
						Eigen::Vector3d normal_vec = dp_dy.cross(dp_dx);
						normal_vec.normalize();

						out_point_cloud.push_back(point_map[y][x]);
						out_normals.push_back(normal_vec);
					}
				}

			}

		}
	}


}

void DepthCameraModel::undistortDepthMap_NoInterpolation(const cv::Mat_<float>& dist_image, cv::Mat_<float>& out_undist_image)
{
	const size_t image_width = dist_image.cols;
	const size_t image_height = dist_image.rows;
	out_undist_image.create(image_height, image_width);

	const Eigen::Matrix3d inv_intrinsics = this->intrinsicMat_.inverse();
	for (size_t undist_y = 0; undist_y < image_height; ++undist_y)
	{
		for (size_t undist_x = 0; undist_x < image_width; ++undist_x)
		{
			Eigen::Vector3d undist_pos_pix(undist_x, undist_y, 1);
			Eigen::Vector3d undist_pos_norm = inv_intrinsics * undist_pos_pix;
			Eigen::Vector3d dist_pos_norm = this->applyRadialDistortion_norm(undist_pos_norm);
			Eigen::Vector3d dist_pos_pix = intrinsicMat_ * dist_pos_norm;

			out_undist_image(undist_y, undist_x) = dist_image((int)(dist_pos_pix[1]), (int)(dist_pos_pix[0]));
		}
	}
}


bool GetNormalFromDepthMap(
	const CameraModel& camera,
	const cv::Mat_<float>& depth_map,
	const double pix_x,
	const double pix_y,
	Eigen::Vector3d& out_normal_vec
)
{
	//const Eigen::Vector3d target_pos_3D 
	//	= camera.computeGlobalPosFromDepth_NoDist(
	//	Eigen::Vector3d(pix_x, pix_y, 1), depth_map((int)pix_y, (int)pix_x));

	//const Eigen::Vector3d next_x_pos_3D 
	//	= camera.computeGlobalPosFromDepth_NoDist(
	//		Eigen::Vector3d(pix_x+1, pix_y, 1), depth_map((int)pix_y, (int)pix_x+1));

	//const Eigen::Vector3d next_y_pos_3D
	//	= camera.computeGlobalPosFromDepth_NoDist(
	//		Eigen::Vector3d(pix_x, pix_y+1, 1), depth_map((int)pix_y+1, (int)pix_x));

	const Eigen::Vector3d target_pos_3D
		= camera.computeGlobalPosFromDepth_NoDist(
			Eigen::Vector3d(pix_x, pix_y, 1), depth_map(floor(pix_y), floor(pix_x)));

	const Eigen::Vector3d next_x_pos_3D
		= camera.computeGlobalPosFromDepth_NoDist(
			Eigen::Vector3d(pix_x + 1, pix_y, 1), depth_map(floor(pix_y), floor(pix_x + 1)));

	const Eigen::Vector3d next_y_pos_3D
		= camera.computeGlobalPosFromDepth_NoDist(
			Eigen::Vector3d(pix_x, pix_y + 1, 1), depth_map(floor(pix_y + 1), floor(pix_x)));


	out_normal_vec = (next_y_pos_3D - target_pos_3D).cross(next_x_pos_3D - target_pos_3D);
	out_normal_vec.normalize();

	return true;
}