#include "../util/util.h"
#include "filters.h"
#include "settings.h"

#include <iostream>

using namespace cv;

template<class T>
class MyLevelFilter : public BaseFilter {
public:
	MyLevelFilter(Mat kernel)
		: kernel(kernel)
	{
		ksize = kernel.size();
		anchor = Point(ksize.width / 2, ksize.height / 2);
	}

	virtual ~MyLevelFilter() { }

	virtual void operator()(const uchar** _src, uchar* _dst, int dststep,
		int dstcount, int width, int cn)
	{
		const T** src = (const T**)_src;
		T* dst = (T*)_dst;

		CV_Assert(cn == 1);
		for (int dstRow = 0; dstRow < dstcount; ++dstRow) {
			for (int dstCol = 0; dstCol < width; ++dstCol) {
				T& accum = dst[dstRow * width + dstCol] = 0;
				T ref = src[dstRow + ksize.height / 2][dstCol + ksize.width / 2];
				for (int i = 0; i < ksize.height; ++i) {
					for (int j = 0; j < ksize.width; ++j) {
						accum += kernel.at<T>(i, j) * fabs(ref - src[dstRow + i][dstCol + j]);
					}
				}
			}
		}
	}

	Mat kernel;
};

Ptr<BaseFilter> getLevelFilter(int type, InputArray filter_kernel)
{
	Mat kernel = filter_kernel.getMat();
	CV_Assert(CV_MAT_TYPE(kernel.type()) == type);
	switch (type)
	{
	case CV_32F:
		return Ptr<BaseFilter>(new MyLevelFilter<float>(kernel));
	case CV_64F:
		return Ptr<BaseFilter>(new MyLevelFilter<double>(kernel));
	default:
		CV_Assert("Unsupported type" == NULL);
		abort();
		break;
	}
}

Ptr<FilterEngine> createLevelFilter(int type, InputArray filter_kernel,
	int _rowBorderType, int _columnBorderType,	const Scalar& _borderValue) 
{
	Ptr<BaseFilter> filter = getLevelFilter(CV_MAT_TYPE(type), filter_kernel);
	return makePtr<FilterEngine>(filter, Ptr<BaseRowFilter>(),
		Ptr<BaseColumnFilter>(), type, type, type,
		_rowBorderType, _columnBorderType, _borderValue);
}

void levelFilter2D(InputArray image, OutputArray dst, int type, InputArray _kernel,
	int _rowBorderType = BORDER_DEFAULT, int _columnBorderType = -1, const Scalar& _borderValue = Scalar())
{
	Mat src = image.getMat();
	src.convertTo(src, type);
	Mat kernel = _kernel.getMat();
	kernel.convertTo(kernel, type);

	Ptr<FilterEngine> filter = createLevelFilter(type, kernel, _rowBorderType, _columnBorderType, _borderValue);

	dst.create(src.size(), type);
	Mat res = dst.getMat();
	filter->apply(src, res);
}

void hessian2D(Mat image, double sigma, OutputArray dxx, OutputArray dxy, OutputArray dyy, bool level) {
	int xmax = int(round(sigma * 3));
	Mat range(1, xmax * 2 + 1, CV_32FC1);
	for (int i = -xmax; i <= xmax; ++i)
		range.at<float>(0, i + xmax) = float(i);
	Mat X, Y, XY;
	repeat(range.reshape(1, 1), range.cols, 1, X);
	repeat(range.reshape(1, 1).t(), 1, range.cols, Y);

	multiply(X, Y, XY);
	multiply(X, X, X);
	multiply(Y, Y, Y);

	Mat weight = (-X - Y) / (2 * sigma * sigma);
	exp(weight, weight);
	Mat DGaussxx, DGaussxy, DGaussyy;
	multiply(1. / (2 * M_PI * pow(sigma, 4)) * (X / (sigma * sigma) - 1), weight, DGaussxx);
	multiply(1. / (2 * M_PI * pow(sigma, 6)) * XY, weight, DGaussxy);
	DGaussyy = DGaussxx.t();

#ifdef _DEBUG_OUT
	Mat tmp;
	cv::hconcat(std::vector<Mat>{ DGaussxx, DGaussxy, DGaussyy }, tmp);
	imwrite("dbg_frangi__2deriv.png", 255 * to01(tmp));
#endif

	if (!level) {
		filter2D(image, dxx, CV_32F, DGaussxx);
		filter2D(image, dxy, CV_32F, DGaussxy);
		filter2D(image, dyy, CV_32F, DGaussyy); 
	}
	else {
		levelFilter2D(image, dxx, CV_32F, DGaussxx);
		levelFilter2D(image, dxy, CV_32F, DGaussxy);
		levelFilter2D(image, dyy, CV_32F, DGaussyy);
	}
#ifdef _DEBUG_OUT
	imwrite("dbg_frangi__im_xy.png", 255 * to01(abs(dxy.getMat())));
	imwrite("dbg_frangi__im_xx.png", 255 * to01(abs(dxx.getMat())));
	imwrite("dbg_frangi__im_yy.png", 255 * to01(abs(dyy.getMat())));
#endif
}

void eig2image(Mat Dxx, Mat Dxy, Mat Dyy, OutputArray l1, OutputArray l2, OutputArray x, OutputArray y) {
	// Compute the eigenvectors of J, v1 and v2
	
	auto diff = Dyy - Dxx;
	Mat tmp;
	sqrt(diff.mul(diff) + 4 * Dxy.mul(Dxy), tmp);
	Mat v2x = 2 * Dxy;
	Mat v2y = diff + tmp;

	// Normalize
	Mat mag;
	sqrt(v2x.mul(v2x) + v2y.mul(v2y), mag); 
	mag.setTo(1, mag == 0);
	v2x = v2x / mag;
	v2y = v2y / mag;

	// The eigenvectors are orthogonal
	Mat v1x = -v2y;
	Mat v1y = v2x;

	// Compute the eigenvalues
	Mat mu1 = 0.5*(Dxx + Dyy + tmp);
	Mat mu2 = 0.5*(Dxx + Dyy - tmp);

	// Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
	Mat check = abs(mu1)>abs(mu2);

	Mat Lambda1; mu1.copyTo(Lambda1); mu2.copyTo(Lambda1, check);
	Mat Lambda2; mu2.copyTo(Lambda2); mu1.copyTo(Lambda2, check);

	Mat Ix; v1x.copyTo(Ix); v2x.copyTo(Ix, check);
	Mat Iy; v1y.copyTo(Iy); v2y.copyTo(Iy, check);
	
	Lambda1.copyTo(l1);
	Lambda2.copyTo(l2);
	Ix.copyTo(x);
	Iy.copyTo(y);

#ifdef _DEBUG_OUT
	REPORT_RANGE(Lambda1);
	REPORT_RANGE(Lambda2);
#endif
}

template<class Op>
void pairwise(Mat src1, Mat src2, OutputArray out, Op op) {
	CV_Assert(src1.size() == src2.size());
	Mat res(src1.size(), CV_32FC1);
	std::transform(src1.begin<float>(), src1.end<float>(), src2.begin<float>(), (float*)res.data, op);
	res.copyTo(out);
}

void frangiFilter2D(Mat src, OutputArray _dst, OutputArray _ori, cv::OutputArray _scale, double sigmaStart, double sigmaStep, double sigmaEnd, double beta, double c, bool level, bool bw) {
	std::vector<Mat> channels;
	split(src, channels);

	beta = 2 * beta * beta;
	c = 2 * c * c;

	Mat res(src.size(), CV_32FC3); // confidence, scale, angle

	Mat Dxx, Dxy, Dyy, Lambda2, Lambda1, Ix, Iy;
	Mat dst(src.size(), CV_32FC1), ori(src.size(), CV_32FC1), scale(src.size(), CV_32FC1);
	dst.setTo(-1);
	ori.setTo(-1);
	scale.setTo(-1);
	for (Mat ch : channels) {
		for (double s = sigmaStart; s <= sigmaEnd; s += sigmaStep) {

			// Make 2D hessian
			hessian2D(ch, s, Dxx, Dxy, Dyy, level);

			// Correct for scale
			Dxx = (s * s) * Dxx;
			Dxy = (s * s) * Dxy;
			Dyy = (s * s) * Dyy;

			// Calculate(abs sorted) eigenvalues and vectors
			eig2image(Dxx, Dxy, Dyy, Lambda1, Lambda2, Ix, Iy);

			// Compute the direction of the minor eigenvector
			Mat angles;
			pairwise(Ix, Iy, angles, [](double x, double y) { return atan2(x, y); });

			// Compute some similarity measures
			Lambda1.setTo(1e-9, Lambda1 == 0);
			Mat Rb = (Lambda1 / Lambda2); 
#ifdef _DEBUG_OUT
			REPORT_RANGE(Rb);
#endif
			Rb = Rb.mul(Rb);
			Mat S2 = Lambda1.mul(Lambda1) + Lambda2.mul(Lambda2);
#ifdef _DEBUG_OUT
			REPORT_RANGE(S2);
#endif

			// Compute the output image
			exp(-Rb / beta, Rb); 
			exp(-S2 / c, S2);
			Mat Ifiltered = Rb.mul(1 - S2);

			// see pp. 45
			if (bw)
				Ifiltered.setTo(0, Lambda2<0);
			else
				Ifiltered.setTo(0, Lambda2>0);
			
			Mat mask = dst < Ifiltered;
			Ifiltered.copyTo(dst, mask);
			angles.copyTo(ori, mask);
			scale.setTo(s, mask);
		}
	}
	if (_dst.needed())
		dst.copyTo(_dst);
	if (_ori.needed())
		ori.copyTo(_ori);
	if (_scale.needed())
		scale.copyTo(_scale);
}
