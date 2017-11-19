#include "volume.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>


void Volume::loadDims() {
	const int* ptr = imageDataRes->GetDimensions();
	for (int i = 0; i < 3; ++i)
		dims[i] = ptr[2-i];
}

Volume::Volume()
	: imageDataRes(vtkSmartPointer<vtkImageData>::New())
{
}

Volume::Volume(vtkImageData* imageData, int dim)
	: imageDataRes(vtkSmartPointer<vtkImageData>::New())
{
	imageDataRes->SetExtent(imageData->GetExtent());
	imageDataRes->SetSpacing(imageData->GetSpacing());
	imageDataRes->SetOrigin(imageData->GetOrigin());
	imageDataRes->AllocateScalars(VTK_FLOAT, dim);

	loadDims();
}

Volume::Volume(const char* fileName)
{
	vtkSmartPointer<vtkXMLImageDataReader> reader =
		vtkSmartPointer<vtkXMLImageDataReader>::New();
	//vtkXMLImageDataReader* reader = vtkXMLImageDataReader::New();
	reader->SetFileName(fileName);
	reader->Update();

	imageDataRes = reader->GetOutput();
	CV_Assert(VTK_FLOAT == imageDataRes->GetScalarType());

	loadDims();
}

float* Volume::getPtr() {
	return (float*)imageDataRes->GetScalarPointer();
}

cv::Mat Volume::getMat() {
	return cv::Mat(3, dims, CV_MAKETYPE(CV_32F, imageDataRes->GetNumberOfScalarComponents()), getPtr());
}

vtkImageData* Volume::getVtkImageData() {
	return imageDataRes.Get();
}

const vtkImageData* Volume::getVtkImageData() const {
	return imageDataRes.Get();
}

void Volume::dumpToVtiFile(const char* fileName) const {
	vtkSmartPointer<vtkXMLImageDataWriter> writer =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();
	writer->SetHeaderType(vtkXMLWriter::UInt64);
	writer->SetFileName(fileName);
	writer->SetInputData(imageDataRes);
	writer->SetCompressorTypeToNone();
	writer->EncodeAppendedDataOff();
	int ret = writer->Write();
	CV_Assert(ret == 1);
}

float Volume::resolution(int dim) const {
	const double* res = imageDataRes->GetSpacing();
	if (dim == -1) {
		CV_Assert(res[0] == res[1] && res[1] == res[2]);
		CV_Assert(res[0] > 0);
		return res[0];
	} 
	CV_Assert(0 <= dim && dim <= 2);
	CV_Assert(res[dim] > 0);
	return res[dim];
}

cv::Vec3f Volume::getOrigin() const {
	return cv::Vec3d(imageDataRes->GetOrigin());
}
