#include "tracker2d.h"

using namespace cv;

struct PointCmp {
	bool operator ()(const Point& lhs, const Point& rhs) const {
		return lhs.x < rhs.x || lhs.x == rhs.x && lhs.y < rhs.y;
	}
};

struct Front {
	Vec2d dir;
	std::vector<Point> Points;
	std::vector<double> Weights;
	std::vector<std::pair<Point, std::vector<int>>> Agg;
};

Front getFront(double gamma, int size)
{
	Point center(size, size);

	int kernelSize = size * 2 + 1;
	Mat kernel(kernelSize, kernelSize, CV_8UC1);
	kernel.setTo(0);
	ellipse(kernel, center, Size(size, size), 0, 270-gamma, 270+gamma, 255, 1, 4);

	imwrite("echo.png", kernel);

	std::map<Point, vector<int>, PointCmp> rel;
	std::vector<Point> points;
	std::vector<double> weights;

	int id = 0;
	for (Point p; p.y < kernelSize; ++p.y){
		for (p.x = 0; p.x < kernelSize; ++p.x){
			if (kernel.at<char>(p)) {
				Point ds = p - center;
				points.push_back(ds);
				weights.push_back(max(0., 1 - fabs(atan2(ds.x, -ds.y)) / (M_PI / 180 * gamma)));
				LineIterator it(kernel, center, p);
				for (int i = 0; i < it .count; ++i, ++it) {
					rel[it.pos() - center].push_back(id);
				}
				++id;
			}
			
		}
	}
	return{ { 0, -1 }, points, weights, { rel.begin(), rel.end() } };
}

Point rotate(Point p, double c, double s) {
	return Point(round(p.x * c + p.y * s), round(-p.x * s + p.y * c));
}

Vec2d rotate(Vec2d p, double c, double s) {
	return Vec2d((p[0] * c + p[1] * s), (-p[0] * s + p[1] * c));
}

Front rotateFront(const Front& f, double angle) {
	double c = cos(angle), s = sin(angle);
	Front ret = f;
	ret.dir = rotate(ret.dir, c, s);
	for (auto& p : ret.Points) {
		p = rotate(p, c, s);
	}
	for (auto& ds : ret.Agg) {
		ds.first = rotate(ds.first, c, s);
	}
	return ret;
}

template<class T>
T getAngleDiff(const cv::Mat& oriMat, size_t oriCnt, Point pos, T alpha) {
	float angle = oriMat.at<ushort>(pos) / (float)oriCnt * M_PI - M_PI_2;
	T diff = T(angle + M_PI * 2) - alpha;
	while (diff > M_PI)
		diff -= T(M_PI); // 0 <= diff < PI
	diff = min(diff, T(M_PI) - diff); // 0 <= diff < PI/2
	T cs = max(T(0), diff - T(M_PI / 180 * 10));
	return cs;
}

double stepHairPoint2D(Point& cur, const Mat& hmap, const Mat& omap, const std::vector<Front>& angleMap, double nu, Vec2d& dir)
{
	Rect box(0, 0, hmap.cols, hmap.rows);
	int sign = dir[0] < 0 ? -1 : 1;
	int alpha = int(acos(-sign * dir[1] / sqrt(1.0 * dir.dot(dir))) / M_PI * angleMap.size() + 0.5);
	if (angleMap.size() == alpha) {
		sign = -1;
		alpha = 0;
	}
	const Front& f = angleMap[alpha];
	std::vector<double> scores(f.Points.size());
	std::vector<int> count(f.Points.size());

	double angle = atan2(dir[1], dir[0]);

	for (auto x : f.Agg) {
		Point next = cur + x.first * sign;
		if (box.contains(next)) {
			float rawHVal = hmap.at<float>(next);
			if (rawHVal >= 0) {
				double h = (rawHVal - nu) / (1 - nu);
				h *= 1 - getAngleDiff(omap, angleMap.size(), next, angle);
				for (auto t : x.second) {
					scores[t] += log(max(0.0, h));
					count[t]++;
				}
			}
		}
	}
	size_t bestId = -1;
	double bestScore = -1;
	for (size_t i = 0; i < scores.size(); ++i) {
		scores[i] = exp(scores[i] / count[i]) * f.Weights[i];
		//scores[i] *= f.Weights[i] / count[i];
		if (box.contains(cur + sign * f.Points[i]) && scores[i] > bestScore) {
			bestId = i;
			bestScore = scores[i];
		}
	}
	if (bestId != -1) {
		cur += sign * f.Points[bestId];
		if (sign*dir[0] * f.Points[bestId].x + sign*dir[1] * f.Points[bestId].y < 0) {
			std::cout << "br" << std::endl;
		}
		dir[0] = sign * f.Points[bestId].x;
		dir[1] = sign * f.Points[bestId].y;
	}
	return bestScore;
}

std::vector<Point> trackHair2D(const Mat& hmap, const Mat& ori, const std::vector<Front>& angleMap, Point start, double th, double nu) {
	std::vector<Point> result;
	Vec2d dir = angleMap[ori.at<unsigned short>(start)].dir;
	for (Point cur = start; stepHairPoint2D(cur, hmap, ori, angleMap, nu, dir) >= th && result.size() < 100;)
		result.push_back(cur);
	std::reverse(result.begin(), result.end());
	result.push_back(start);
	if (result.size() > 1) {
		Point diff = result[result.size() - 1] - result[result.size() - 2];
		dir[0] = diff.x;
		dir[1] = diff.y;
	} else {
		dir = -angleMap[ori.at<unsigned short>(start)].dir;
	}
	for (Point cur = start; stepHairPoint2D(cur, hmap, ori, angleMap, nu, dir) >= th && result.size() < 100;)
		result.push_back(cur);
	return result;
}

double avg(const Mat& mask, const std::vector<Point>& track) {
	double sum = 0;
	int count = 0;
	for (size_t i = 1; i < track.size(); ++i) {
		LineIterator iter(mask, track[i - 1], track[i]);
		count += iter.count;
		for (int j = 0; j < iter.count; ++j, ++iter)
			if (mask.at<char>(iter.pos()))
				++sum;
	}
	return sum / count;
}

//int main(int argc, const char* argv[]) {
// command line argument:
// "C:\Users\shangxuanu\Desktop\dmitry_hair-master\dmitry_hair-master\hair\detector\330011.yml" "C:\Users\shangxuanu\Desktop\dmitry_hair-master\dmitry_hair-master\hair\detector\330011_ori.png" res.txt_2dlines --non-max 

bool tracker2d(std::string input_yml, std::string ori_png, std::string output_fn) {
	bool non_max = true;

	std::string hmapPath = input_yml;
	std::string oriPath = ori_png;
	std::string outPath = output_fn;
	
	int step = 10;
	double gamma = 60;
	double th = 0.1, nu = 0, seedThreshold = 0.2;
	int lambda = 3, minLen = 20;
	std::vector<int> dbgId;
	
	/*namespace po = boost::program_options;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("i", po::value(&hmapPath)->required(), "hmaps or gabor result")
		("o", po::value(&oriPath)->required(), "orientation map result")
		("r", po::value(&outPath)->required(), "result path")
		("s", po::value(&step)->default_value(10), "step size")
		("g", po::value(&gamma)->default_value(60), "gamma")
		("t", po::value(&th)->default_value(0.1), "threshold")
		("st", po::value(&seedThreshold)->default_value(0.2), "seed threshold")
		("nu", po::value(&nu)->default_value(0), "nu")
		("la", po::value(&lambda)->default_value(3), "hair thikness, px")
		("debug", po::value(&dbgId), "debug for hair ids")
		("min-hair-length", po::value(&minLen)->default_value(20), "min hair length, px")
		("non-max", "non-max for seeds")
		("2d", "set z = 0")
		;

	po::positional_options_description p;
	p.add("i", 1);
	p.add("o", 1);
	p.add("r", 1);

	po::variables_map vm;

	try {
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl
			<< desc << std::endl;
		return -1;
	}*/



	try {

		Mat hmap = unifiedRead(hmapPath);
		if (!hmap.data) {
			fprintf(stderr, "Cannot read file %s\n", hmapPath);
			return -1;
		}
		if (hmap.type() == CV_16UC1)
			hmap.convertTo(hmap, CV_32F, 1.0 / USHRT_MAX);
		else
			CV_Assert(hmap.type() == CV_32FC1);

		Mat ori = imread(oriPath, CV_LOAD_IMAGE_UNCHANGED);
		if (!ori.data) {
			fprintf(stderr, "Cannot read file %s\n", oriPath);
			return -1;
		}

		if (ori.type() != CV_16U) {
			fprintf(stderr, "Expected type is uint16");
			return -1;
		}

		std::vector<unsigned short> angles(ori.begin<unsigned short>(), ori.end<unsigned short>());
		std::sort(angles.begin(), angles.end());
		angles.erase(std::unique(angles.begin(), angles.end()), angles.end());
		ori = ori * angles.size() / USHRT_MAX;
		Front f = getFront(gamma, step);
		std::vector<Front> angleMap(angles.size());
		//std::vector<Front> angleMap(angles.back() + 1);
		for (int i = 0; i < angles.size(); ++i)
			angleMap[i] = rotateFront(f, -angles[i] / 65536. * M_PI);
		//angleMap[angles[i]] = rotateFront(f, -angles[i] / 65536. * M_PI);

		Mat seedMat(hmap.size(), CV_32FC1);
		seedMat = 0;
		hmap.copyTo(seedMat, hmap >= seedThreshold);
		if (non_max) {
			Mat oriTmp;
			ori.convertTo(oriTmp, CV_32SC1);
			nonMax(oriTmp, seedMat, lambda, angles.size(), seedMat, 1);
		}
		imwrite("seeds.png", to01(seedMat) * 255);

		std::vector<std::pair<float, Point>> seeds;
		for (int y = 0; y < hmap.rows; ++y) {
			for (int x = 0; x < hmap.cols; ++x) {
				if (seedMat.at<float>(y, x) >= seedThreshold) {
					seeds.emplace_back(hmap.at<float>(y, x), Point(x, y));
				}
			}
		}
		struct cmp {
			bool operator()(const std::pair<float, Point>& lhs, std::pair<float, Point>& rhs) {
				return lhs.first > rhs.first;
			}
		};
		std::sort(seeds.begin(), seeds.end(), cmp());

		std::vector<std::pair<std::vector<Point>, size_t>> tracks;
		for (size_t id = 0; id < seeds.size(); ++id) {
			if (!dbgId.empty() && std::find(dbgId.begin(), dbgId.end(), id) == dbgId.end())
				continue;
			auto& pair = seeds[id];
			if (hmap.at<float>(pair.second) >= 0) {
				auto& seed = pair.second;
				std::vector<Point> hair = trackHair2D(hmap, ori, angleMap, seed, th, nu);
				polylines(hmap, hair, false, -1, lambda);
				if (!dbgId.empty()) {
					circle(hmap, seed, 3, -2);
					imwrite("hmap.png", to01(hmap) * 255);
				}
				if (hair.size() * step >= minLen + step)
					tracks.emplace_back(hair, id);
			}
		}
		std::ofstream ofs(outPath);

		Mat hairs(hmap.size(), CV_8UC1);
		hairs.setTo(0);
		Mat hairsIdx(hmap.size(), CV_16UC1);
		hairsIdx.setTo(0);
		int written = 0;
		for (size_t i = 0; i < tracks.size(); ++i) {
			if (tracks[i].first.size() > 1) {
				polylines(hairsIdx, tracks[i].first, false, i, 1);
				polylines(hairs, tracks[i].first, false, 255, 1);
				dumpTrack(ofs, tracks[i].first, tracks[i].second, seeds[tracks[i].second].first) << std::endl;
				++written;
			}
		}
		std::cout << written << " tracks out of " << tracks.size() << " saved" << std::endl;
		imwrite("hmap.png", to01(hmap) * 255);
		imwrite("hairs.png", hairs);
		imwrite("hairs_idx.png", hairsIdx);
		return true;
	}
	catch (const std::exception& ex) {
		std::cerr << hmapPath << ": " << ex.what() << std::endl;
		return false;
	}
}
