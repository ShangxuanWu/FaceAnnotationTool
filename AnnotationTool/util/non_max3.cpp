#include "util.h"
#include "convolution3.h"

using namespace cv;

static const float s2 = sqrt(1 / 2.0f);
static const float s3 = sqrt(1 / 3.0f);

static const int MAJOR_DIR_NUM = 13;
static const Vec3f major_dirs[MAJOR_DIR_NUM] =
{
	// directions along the axis, there are 3 of them
	Vec3f(1, 0, 0), Vec3f(0, 1, 0), Vec3f(0, 0, 1),
	// directions that lie within the 3D plane of two axis, there are 6 of them
	Vec3f(s2, s2, 0), Vec3f(0, s2, s2), Vec3f(s2, 0, s2),
	Vec3f(s2, -s2, 0), Vec3f(0, s2, -s2), Vec3f(s2, 0, -s2),
	// directions that are equally in between three axis, there are 4 of them
	Vec3f(s3, s3, s3), Vec3f(s3, s3, -s3), Vec3f(s3, -s3, s3), Vec3f(s3, -s3, -s3)
};

static const int NUM_DIR_3D = 26;
static Vec3i neighbor3d[NUM_DIR_3D] =
{
	// directions along the axis, there are 6 of them
	Vec3i(0, 0, 1), Vec3i(0, 1, 0), Vec3i(1, 0, 0),
	Vec3i(0, 0, -1), Vec3i(0, -1, 0), Vec3i(-1, 0, 0),
	// directions that lie within the 3D plane of two axis, there are 12 of them
	Vec3i(0, 1, 1), Vec3i(0, 1, -1), Vec3i(0, -1, 1), Vec3i(0, -1, -1),
	Vec3i(1, 0, 1), Vec3i(1, 0, -1), Vec3i(-1, 0, 1), Vec3i(-1, 0, -1),
	Vec3i(1, 1, 0), Vec3i(1, -1, 0), Vec3i(-1, 1, 0), Vec3i(-1, -1, 0),
	// directions that are equally in between three axis, there are 8 of them
	Vec3i(1, 1, 1), Vec3i(1, 1, -1), Vec3i(1, -1, 1), Vec3i(1, -1, -1),
	Vec3i(-1, 1, 1), Vec3i(-1, 1, -1), Vec3i(-1, -1, 1), Vec3i(-1, -1, -1),
};

vector<vector<Vec3i>> initCrossSection() {
	vector<vector<Vec3i>> cross_section(MAJOR_DIR_NUM);
	// Setting offsets that are perpendicular to dirs
	for (int i = 0; i < MAJOR_DIR_NUM; i++)
	{
		for (int j = 0; j < NUM_DIR_3D; j++)
		{
			// multiply the two directions
			float temp = major_dirs[i].dot(neighbor3d[j]);
			// the temp is 0, store this direciton
			if (fabs(temp) < 1.0e-5)
			{
				cross_section[i].push_back(neighbor3d[j]);
			}
		}
	}
	return cross_section;
}

// cross section for the major orientation
static vector<vector<Vec3i>> cross_section = initCrossSection();

#define _DEBUG_OUT

void nonMax3(const Mat& conf, const Mat& dir, double sigma, OutputArray dst) {
	CV_Assert(conf.dims == 3);

	int dims[] = { conf.size[0], conf.size[1], conf.size[2] };

	dst.create(3, dims, CV_32FC1);
	Mat result = dst.getMat();

	Mat data;
	if (sigma)
		gaussianBlur3(conf, sigma, size_t(round(sigma)) * 6 + 1, data);
	else
		data = conf;

	
#pragma omp parallel for
	for (int z = 0; z < data.size[0]; ++z) {
		for (int y = 0; y < data.size[1]; ++y) {
			for (int x = 0; x < data.size[2]; ++x) {
				bool isMax = true;
				const Vec3f& cur_dir = dir.at<Vec3f>(z, y, x);
				
				float len = cur_dir.dot(cur_dir);
				if (len == 0)
					continue;

				CV_Assert(fabs(len - 1) < 1e-5);

				float ref = data.at<float>(z, y, x);

				int mdi = 0; // major direction id
				float max_dot_product = 0;
				for (int di = 0; di<MAJOR_DIR_NUM; di++)
				{
					float current_dot_product = fabs(cur_dir.dot(major_dirs[di]));
					if (max_dot_product < current_dot_product)
					{
						max_dot_product = current_dot_product;
						mdi = di;// update the major direction id
					}
				}

				for (size_t i = 0; i<cross_section[mdi].size(); i++)
				{
					size_t ox = x + cross_section[mdi][i][0];
					size_t oy = y + cross_section[mdi][i][1];
					size_t oz = z + cross_section[mdi][i][2];
					if (ox < data.size[2] && oy < data.size[1] && oz <  data.size[0] 
						&& ref < data.at<float>(oz, oy, ox))
					{
						isMax = false;
						break;
					}
				}

				if (isMax)
				{
					result.at<float>(z, y, x) = conf.at<float>(z, y, x);
				}
			}
		}
	}

}
