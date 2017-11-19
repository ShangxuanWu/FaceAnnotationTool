#pragma once

#include <opencv2/core.hpp>

#include <vtkSmartPointer.h>

class vtkImageData;

// data stored in order z, y, x. Fastest changing index is x
class Volume {
public:
	Volume();

	explicit Volume(vtkImageData* imageData, int dim = 1);

	explicit Volume(const char* fileName);

	const int* getDims() const { return dims; }

	float* getPtr();

	cv::Mat getMat();

	vtkImageData* getVtkImageData();

	const vtkImageData* getVtkImageData() const;

	void dumpToVtiFile(const char* fileName) const;

	float resolution(int dim = -1) const;

	cv::Vec3f getOrigin() const;

	cv::Vec3f fromMatCoo(const cv::Vec3f& coo) const {
		cv::Vec3f result(coo);
		std::swap(result(0), result(2)); // volume coordinate order is z, y, x
		return result * resolution() + getOrigin();
	}

	cv::Vec3i toMatCoo(const cv::Vec3f& coo) const {
		cv::Vec3f result = (coo - getOrigin()) / resolution();
		std::swap(result(0), result(2));
		return result;
	}

	template<class Vec3>
	cv::Vec3i toMatCoo(const Vec3& coo) const {
		return toMatCoo({ (float)coo[0], (float)coo[1], (float)coo[2]});
	}

private:
	vtkSmartPointer<vtkImageData> imageDataRes;
	int dims[3];

	Volume(const Volume&);

	void loadDims();
};