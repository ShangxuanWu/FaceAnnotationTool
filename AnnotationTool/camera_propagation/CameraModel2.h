#ifndef _CAMERA_MODEL2_H_
#define _CAMERA_MODEL2_H_

#include <string>

//#include "CalibrationIO.h"
//#include "OpenCV_util.h"

#include <opencv2/opencv.hpp>

//#include "EigenUtil.h"
#include <Eigen/Eigen>

enum CameraType
{
	Camera_Undefined = -1,
	Camera_RGB = 0,
	Camera_IRDEPTH,
	NumCameraTypes
};

class CameraModel
{
public:
//	static string CalibParamNames;

	CameraModel();
//	CameraModel(const std::string& calibration_filepath);
	CameraModel(const CameraModel& copy);
	virtual ~CameraModel(){};

	//bool loadCalibrationFile(const std::string& calib_filepath);

	virtual bool loadCalibrationFile(const std::string& calib_filepath);
	virtual bool writeCalibrationFile(const std::string& calib_filepath) const;

	// extrinsics
	void setExtrinsicMat(const Eigen::MatrixXd& extrinsic_mat);
	void setExtrinsicParams(const Eigen::Matrix3d& global_orient, const Eigen::Vector3d& global_pos);
	Eigen::MatrixXd extrinsicMat() const;
	cv::Mat extrinsicMat_cv() const;
	const Eigen::Matrix3d& globalOrient_mat() const { return globalOrient_mat_; }
	const Eigen::Vector3d& globalPos() const { return globalPos_; }
	Eigen::Vector3d t() const { return -1. * globalOrient_mat_ * globalPos_; }

	virtual bool loadImage(const std::string& filepath){return false;}

	// intrinsics
	Eigen::Matrix3d intrinsicMat() const { return intrinsicMat_; }
	void setIntrinsicParams(double focal_x, double focal_y, double pp_x, double pp_y);
	void setIntrinsicMat(const Eigen::Matrix3d& intrinsic_mat);// { intrinsicMat_ = intrinsic_mat; }
	double focalLengthX() const { return intrinsicMat_(0, 0); }
	double focalLengthY() const { return intrinsicMat_(1, 1); }
	double principalPointX() const { return intrinsicMat_(0, 2); }
	double principalPointY() const { return intrinsicMat_(1, 2); }

	cv::Mat& intrinsicMat_cv() { return intrinsicMat_cv_; }
	const cv::Mat& intrinsicMat_cv() const { return intrinsicMat_cv_; }

	// lens distortion parameters
	void setDistortionParams(double params[5])
	{
		distortionParams_.resize(5);
		distortionParams_cv_.create(5, 1, CV_64FC1);

		for (int i = 0; i < 5; ++i)
		{
			distortionParams_[i] = params[i];
			distortionParams_cv_.at<double>(i, 0) = params[i];
		}
	}
	void setDistortionParams(const Eigen::VectorXd& dist_params) 
	{ 
		distortionParams_ = dist_params; 

		distortionParams_cv_.create(5, 1, CV_64FC1);

		for (int i = 0; i < 5; ++i)
		{
			distortionParams_cv_.at<double>(i, 0) = dist_params[i];
		}

	}
	void setDistortionParams(double k1, double k2, double p1, double p2, double k3)
	{
		Eigen::VectorXd dist_params( 5 );
		dist_params[0] = k1;
		dist_params[1] = k2;
		dist_params[2] = p1;
		dist_params[3] = p2;
		dist_params[4] = k3;

		this->setDistortionParams(dist_params);

	}
	Eigen::VectorXd& distortionParams() { return distortionParams_; }
	const Eigen::VectorXd& distortionParams() const { return distortionParams_; }
	cv::Mat& distortionParams_cv() { return distortionParams_cv_; }
	const cv::Mat& distortionParams_cv() const { return distortionParams_cv_; }

	void setImageResolution(size_t image_width, size_t image_height);
	size_t imageWidth() const { return imageWidth_; }
	size_t imageHeight() const { return imageHeight_; }

	int cameraType() const { return cameraType_; }

	Eigen::Vector3d globalRayDirOfPixel(const Eigen::Vector2d& pixel_pos) const;
	Eigen::Vector3d globalRayDirOfPixel_NoDist(const Eigen::Vector2d& pixel_pos) const;
	Eigen::Vector3d backProjectTo3DWithDistance(const Eigen::Vector2d& pixel_pos, const double distance) const;

	Eigen::Vector3d computeGlobalPosFromDepth_NoDist(const Eigen::Vector2d& pixel_pos, const double depth) const;
	Eigen::Vector3d computeGlobalPosFromDepth_NoDist(const Eigen::Vector3d& pixel_pos_homo, const double depth) const;
	double computeDepth(const Eigen::Vector3d& global_pos_3D) const;

	bool isPointInFrontOfCamera(const Eigen::Vector3d& pos_3D, const double thresh_local_z = 0) const;

	Eigen::Vector3d convertGlobalPosToCameraLocalPos(const Eigen::Vector3d& global_pos_3D) const;

	Eigen::Vector3d projectPoint3DToImageNoDist_pix(const Eigen::Vector3d& pos_3D) const;
	Eigen::Vector3d projectPoint3DToImageWithDist_pix(const Eigen::Vector3d& pos_3D) const;

	Eigen::Vector3d applyRadialDistortion_norm(const Eigen::Vector3d& undist_pos_norm) const;
	
	Eigen::Vector3d undistortPixelPos_norm(const Eigen::Vector3d& dist_pos_norm) const;
	Eigen::Vector3d undistortPixelPos_pix(const Eigen::Vector3d& dist_pos_pix) const;

	void setName(const std::string& camera_name) { cameraName_ = camera_name; }
	const std::string name() const { return cameraName_; }
	const int nameAsInteger() const { return atoi(cameraName_.c_str()); }

protected:
	Eigen::Matrix3d intrinsicMat_;
	Eigen::VectorXd distortionParams_;

	// for undistortion using opencv
	cv::Mat intrinsicMat_cv_; // for undistortion
	cv::Mat distortionParams_cv_;

	// rotation in global coordinate system
	Eigen::Matrix3d globalOrient_mat_;
	// position in global coordinate system
	Eigen::Vector3d globalPos_;


	std::string cameraName_;
	int cameraType_;
	int cameraID_;
	size_t imageWidth_;
	size_t imageHeight_;
};

class ColorCameraModel : public CameraModel
{
public:
	ColorCameraModel();
	ColorCameraModel(const ColorCameraModel& copy);
	~ColorCameraModel();

	bool loadImage(const std::string& image_filepath);
	bool writeImage(const std::string& image_filepath) const;
	//bool applyUndistortion();
	cv::Mat colorImage() const { return colorImage_; }
	void setColorImage(const cv::Mat& color_image) { colorImage_ = color_image; }

//	bool loadImageSeqList(const std::string& directory, const std::string& image_list_filename);

//	size_t numImages() const { return imageFilePathSeq_.size(); }

private:
	cv::Mat colorImage_;

//	std::vector<std::string> imageFilePathSeq_;

};

class DepthCameraModel : public CameraModel
{
public:
	DepthCameraModel();
//	DepthCameraModel(const std::string& calibration_filepath);
	DepthCameraModel(const DepthCameraModel& copy);
	~DepthCameraModel(){};

	bool loadImage(const std::string& png_filepath);
	bool writePseudoDepthImage(const std::string& png_filepath, const double scale ) const;
	//bool applyUndistortion();
	cv::Mat_<float> depthImage_m() const { return depthImage_m_; }

	//void applySmoothing();
	//void computeNormalMap();

	//void convertDepthMapToPointCloud(cv::Mat_<double>& point_pos); //std::vector<cv::Point3f>& point_pos);
	//cv::Mat_<double> normalMap() const { return normalMap_; }

	//void setValidDepthRange(double min_depth_m, double max_depth_m)
	//{
	//	minDepth_m_ = min_depth_m;
	//	maxDepth_m_ = max_depth_m;
	//}

	//bool hasValidDepthAt(const size_t x, const size_t y) const { return validDepthMask_(y, x); }

	//cv::Mat_<double> computePointPosAt(size_t pix_x, size_t pix_y) const;

	//bool loadDepthImageSeqList(const std::string& directory, const std::string& image_list_filename);
	//bool loadIRImageSeqList(const std::string& directory, const std::string& image_list_filename);

	void reconstructPointCloud(const cv::Mat_<float>& depth_map_m, const double max_depth_m, std::vector<Eigen::Vector3d>& out_point_cloud);
	void reconstructPointCloudNormals(
		const cv::Mat_<float>& depth_map_m, const double max_depth_m,
		const int boundary_margin_pix,
		std::vector<Eigen::Vector3d>& out_point_cloud,
		std::vector<Eigen::Vector3d>& out_normals);


	void undistortDepthMap_NoInterpolation(const cv::Mat_<float>& dist_image, cv::Mat_<float>& out_undist_image);

private:
	cv::Mat_<float> depthImage_m_;
//	cv::Mat_<cv::Vec3d> normalMap_;

	inline size_t array1DIdxFromPixXY_(size_t x, size_t y)
	{
		return (y * imageWidth_ + x);
	}

	//cv::Mat_<bool> validDepthMask_;
	double minDepth_m_;
	double maxDepth_m_;
	inline bool hasDepthInRange_(double target_depth)
	{
		return (minDepth_m_ < target_depth && target_depth < maxDepth_m_);
	}
	//void createValidDepthMask_();

	std::vector<std::string> depthImageFilePathSeq_;
	std::vector<std::string> irImageFilePathSeq_;

};

bool GetNormalFromDepthMap(
	const CameraModel& camera,
	const cv::Mat_<float>& depth_map,
	const double pix_x,
	const double pix_y,
	Eigen::Vector3d& out_normal_vec
);

#endif  // _CAMERA_MODEL2_H_