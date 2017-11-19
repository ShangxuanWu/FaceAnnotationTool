#include "detector.h"

using namespace cv;

void frangiForHair(Mat src, int directions, double sigma, double lambda, OutputArray respOut, OutputArray angleOut) {
	int kernel_size = int(lambda * 3) | 1;

	cv::cvtColor(src, src, CV_BGR2HSV);
	src.convertTo(src, CV_32F);

	Mat ch[3];

	split(src, ch);
	imwrite("h.png", ch[0]);
	imwrite("s.png", ch[1]);
	imwrite("v.png", ch[2]);
	Mat tmp, scale;
	char buf[1024];
	for (int i = 2; i < 3; ++i) {
		frangiFilter2D(ch[i], tmp, noArray(), scale, 1.5, 0.5, 1.6, 0.5, 12);
		sprintf(buf, "frangi_%d.png", i);
		imwrite(buf, 255 * to01(tmp));
		sprintf(buf, "scale_%d.png", i);
		imwrite(buf, 255 * to01(scale));
	}
}

int frangiBased(Mat src, Mat mask, OutputArray confidence, OutputArray hairMap, OutputArray orientation, double sigma, double lambda) {
	Mat dst, mx, angle;

	int imax = 18;
	double upper = 0.5;
	double lower = 0.2;

	frangiForHair(src, imax, sigma, lambda, dst, angle);

	return -1;
}

void gaborForHair(Mat src, Mat mask, int directions, double sigma, double lambda, OutputArray respOut, OutputArray angleOut) {
	int kernel_size = int(lambda * 3) | 1;

	cv::cvtColor(src, src, CV_BGR2HSV);
	src.convertTo(src, CV_32F);

	Mat ch[3];

#ifdef _DEBUG_OUT
	split(src, ch);
	imwrite("s.png", ch[1]);
	imwrite("v.png", ch[2]);
#endif

	// crop by mask
	Rect rect;
	Size osize = src.size();
	if (!mask.empty()) {
		cv::Mat nnz;
		findNonZero(mask, nnz);
		if (nnz.empty()) {
			respOut.create(osize, CV_32FC1);
			Mat respMat = respOut.getMat();
			respMat = 0.0;
			angleOut.create(osize, CV_32SC1);
			Mat angleMat = angleOut.getMat();
			angleMat = 0.0;
			return;
		}
		rect = boundingRect(nnz);
		rect.x -= kernel_size;
		rect.y -= kernel_size;
		rect.width += kernel_size * 2;
		rect.height += kernel_size * 2;
		rect = rect & Rect({}, src.size());
		src = src(rect);
	}

	Mat dst = Mat::zeros(src.size(), CV_32FC1);
	Mat angle = Mat::zeros(src.size(), CV_32SC1);
	for (int i = 0; i < directions; ++i) {
		double ang = i *  M_PI / directions;

		//printf(".");

		/// Update kernel size for a normalized box filter
		Mat kernel = getGaborKernel(Size(kernel_size, kernel_size), sigma, ang, lambda, 1, 0);
		kernel = kernel - 1. * sum(kernel)[0] / kernel.total();

		//printf("%lf\n", sum(abs(kernel))[0]);

		/// Apply filter
		Mat cur;
		filter2D(src, cur, -1, kernel);
		split(cur, ch);
		//cur = abs(ch[1]) + abs(ch[2]);

		cur = -ch[2];

		//Mat curAngle(angle.size(), angle.type());
		//curAngle.setTo(i);
		//nonMax(curAngle, cur, lambda, directions, cur, 1);

		// max
		Mat msk = cur > dst;

		cur.copyTo(dst, msk);
		angle.setTo(i, msk);
	}

	if (!mask.empty()) {
		respOut.create(osize, CV_32FC1);
		angleOut.create(osize, CV_32SC1);
		Mat respMat = respOut.getMat();
		respMat = 0.0;
		Mat angleMat = angleOut.getMat();
		angleMat = 0.0;
		dst.copyTo(respMat(rect));
		angle.copyTo(angleMat(rect));
	}
	else {
		dst.copyTo(respOut);
		angle.copyTo(angleOut);
	}
}

int gaborBased(Mat src, Mat mask, OutputArray confidence, OutputArray hairMap, OutputArray orientation, double sigma, double lambda, double upper, double lower) {
	/// Declare variables
	Mat dst, mx, angle;

	int imax = GABOR_NUMBER_OF_DIRECTION;

	gaborForHair(src, mask, imax, sigma, lambda, dst, angle);
	
#ifdef _DEBUG_OUT // added by Shangxuan, but not working
	std::cout << "print dst" << std::endl;
	namedWindow("dst", WINDOW_AUTOSIZE);// Create a window for display.
	imshow("dst", dst);
	waitKey(0);
#endif

	if (!mask.empty()) {
		dst.setTo(0, mask == 0);
	}

#ifdef _DEBUG_OUT
	Mat tmp;
	to01(dst).convertTo(tmp, CV_16U, 65535);
	imwrite("gabor_response.png", tmp);
	imwrite("gabor_orientation.png", angle2rgb(angle, 180. / imax, dst));
#endif
	nonMax(angle, dst, 0, imax, dst);

	Mat bin = dst / 256;

#ifdef _DEBUG_OUT
	to01(dst).convertTo(tmp, CV_16U, 65535);
	imwrite("gabor_response_nm.png", tmp);

	imwrite("bin_mask_lower.png", bin >= lower);
	imwrite("bin_mask_upper.png", bin >= upper);
#endif 

	cannyThreshold(bin, lower, upper, bin);
	
#ifdef _DEBUG_OUT
	imwrite("bin_mask.png", bin);

	dst = to01(dst);
	imwrite("y.png", angle2rgb(angle, 180. / 18, dst > 0.07));
#endif

	//dilate(bin, bin, Mat::ones(2, 2, CV_8UC1));

	Mat dist;
	distanceTransform(bin == 0, dist, CV_DIST_L2, CV_DIST_MASK_PRECISE);

	Mat hair = 1 / (1 + dist);

	hair.copyTo(hairMap);

	dst = to01(dst)*255;
	dst.copyTo(confidence);

	angle.convertTo(angle, CV_16U, 65536. / imax);
	angle.copyTo(orientation);


	return 0;
}

void printHelp(const char* file) {
	fprintf(stderr, "Usage: %s <infile> <mask> <out conf> <out angle> gabor <param1> <param2>... \n", file);
}

#include <time.h>
#include <boost/program_options.hpp>


// int main(int argc, char** argv)
// how they call detector main in command line:
// "C:\Users\shangxuanu\Downloads\image0000.png" none 330011.yml 330011_ori.png gabor_only 2 6 --inpaint-over-exposed 130
bool detector(std::string original_image_path, std::string mask_path, std::string output_yml_name, std::string output_ori_png_name) {
	
	/*namespace po = boost::program_options;
	po::options_description desc("Allowed options");*/

	std::string srcPath = original_image_path,
		maskPath = mask_path, 
		outputScorePath = output_yml_name, 
		outputOriPath = output_ori_png_name, 
		type = "gabor_only";
	
	double sigma = 2, // sigma = 0 means no sigma
		lambda = 6, // lambda = 0 means no lambda
		upper = GABOR_CANNY_UPPER,
		lower = GABOR_CANNY_LOWER;
	bool inpaint_over_exposed = true;
	double overexposed = 130;
	
	/*
	desc.add_options()
		("input", po::value(&srcPath)->required(), "input image")
		("mask", po::value(&maskPath)->required(), "mask")
		("out-score", po::value(&outputScorePath)->required(), "output Score Path")
		("out-ori", po::value(&outputOriPath)->required(), "output Ori Path")
		("type", po::value(&type)->required(), "type")
		("sigma", po::value(&sigma), "sigma")
		("lambda", po::value(&lambda), "lambda")
		("upper", po::value(&upper)->default_value(GABOR_CANNY_UPPER), "Canny upper threshold")
		("lower", po::value(&lower)->default_value(GABOR_CANNY_LOWER), "Canny lower threshold")
		("inpaint-over-exposed", po::value(&overexposed), "removes over exposed region (pixels with intensity >= ARG)")
		;

	po::positional_options_description p;

	p.add("input", 1);
	p.add("mask", 1);
	p.add("out-score", 1);
	p.add("out-ori", 1);
	p.add("type", 1);
	p.add("sigma", 1);
	p.add("lambda", 1);

	po::variables_map vm;
	*/
	bool testRun = true;

	try {
		/*po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);

		po::notify(vm);*/

		AutoTimer timer(srcPath.c_str());

		Mat src = imread(srcPath);

		Mat confidence, hairmap, orientation;

		Mat mask = imread(maskPath, CV_LOAD_IMAGE_UNCHANGED);
		
		Mat ovrx;
		if (inpaint_over_exposed) {
			Mat grey, m;
			cvtColor(src, grey, CV_BGR2GRAY);
			ovrx = grey > overexposed;
			dilate(ovrx, m, getStructuringElement(MORPH_ELLIPSE, { 5, 5 }));
			if (!mask.empty())
				m &= mask;
			inpaint(src, m, src, 3, INPAINT_NS);
		}

		if (!src.data)
		{
			fprintf(stderr, "Empty input!\n");
			printHelp(original_image_path.c_str());
			return -1;
		}

		CV_Assert(mask.empty() || mask.elemSize() == 1);

		if (mask.empty())
			std::cerr << "No mask used" << std::endl;
		else
			erode(mask, mask, getStructuringElement(MORPH_ELLIPSE, { 5, 5 }));

		if (type == std::string("gabor")) {
			if ( sigma == 0 || lambda == 0) {
				throw std::runtime_error("gabor_only: needs sigma and lambda parameters");
			}

			int ret = gaborBased(src, mask, confidence, hairmap, orientation, sigma, lambda, upper, lower);

			if (ret)
				return ret;

			hairmap.convertTo(hairmap, CV_16U, 65535);
			imwrite(outputScorePath, hairmap);
			imwrite(outputOriPath, orientation);

		}
		else if (type == std::string("frangi")) {
			if ( sigma == 0 || lambda == 0) {
				throw std::runtime_error("gabor_only: needs sigma and lambda parameters");
			}

			int ret = frangiBased(src, mask, confidence, hairmap, orientation, sigma, lambda);

			if (ret)
				return ret;

			hairmap.convertTo(hairmap, CV_16U, 65535);
			imwrite(outputScorePath, hairmap);
			imwrite(outputOriPath, orientation);

		}
		else if (type == std::string("gabor_only")) {
			if (sigma == 0 || lambda == 0) {
				throw std::runtime_error("gabor_only: needs sigma and lambda parameters");
			}
			int imax = GABOR_NUMBER_OF_DIRECTION;

			Mat dst;
			Mat angle;
			gaborForHair(src, mask, imax, sigma, lambda, dst, angle);

			dst /= 255;

			if (!ovrx.empty())
				dst.setTo(0, ovrx);

			if (mask.empty())
				floatWrite(outputScorePath, dst);
			else {
				dst.setTo(0, mask == 0);
				floatWrite(outputScorePath, dst, mask);
			}

			angle.convertTo(angle, CV_16U, 65536. / imax);
			imwrite(outputOriPath, angle);
		}
		else {
			throw std::runtime_error("unknown type");
		}
	}
	/*catch (const po::error& ex) {
		std::cerr << "Error :" << ex.what() << "\n";
		std::cerr << desc << "\n";
		return -1;
	}*/
	catch (const std::exception& ex) {
		std::cerr << "Error :" << ex.what() << "\n";
		return -1;
	}
}
