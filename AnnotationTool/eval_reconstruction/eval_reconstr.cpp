#include "../util/util.h"
#include "../util/geom.h"
#include "../util/cutil.h"
#include "../cw_lib/CCameraArray.h"


#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "successive_shortest_path_nonnegative_weights.hpp"

#include <string>
#include <fstream>

namespace po = boost::program_options;

using namespace cv;

std::pair<Mat, Point> computeBand(const std::vector<CVec2f>& _pred, float tol) {
	std::vector<Point> pred;
	for (auto& p : _pred)
		pred.push_back(vector_cast<Point>(p));

	Rect rect = boundingRect(pred);
	cv::Mat mask(rect.size(), CV_8UC1);
	mask = 0;
	for (auto& p : pred)
		p -= rect.tl();
	polylines(mask, pred, false, 255, max(1.f, tol));
	return{ mask, rect.tl() };
}

float computeScore(const std::vector<CVec2f>& _gt, const std::pair<Mat, Point>& mask) {
	std::vector<Point> gt;
	for (auto& p : _gt)
		gt.push_back(vector_cast<Point>(p));
	for (auto& p : gt)
		p -= mask.second;

	int count = 0;
	double total = 1;
	Point prev(-1, -1);
	for (size_t i = 1; i < gt.size(); ++i) {
		LineIterator it(mask.first, gt[i - 1], gt[i]);
		for (int k = 0; k < it.count; ++k, ++it) {
			if (**it && prev != it.pos())
				++count;
			prev = it.pos();
		}
		Point diff = gt[i - 1] - gt[i];
		int curTotal = max(abs(diff.x), abs(diff.y));
		CV_Assert(it.count <= curTotal + 1);
		total += curTotal;
	}
	return count / (float)total;
}

void evaluate(const std::vector<std::vector<CVec2f>>& detected, const std::vector<std::vector<CVec2f>>& gt, float proximityTolerance) {
	using namespace boost;

	typedef adjacency_list<vecS, vecS, directedS> BaseGraph;

	typedef adjacency_list<vecS, vecS, directedS, no_property,
		property<edge_weight_t, float,
		property<edge_capacity_t, int,
		property<edge_residual_capacity_t, int,
		property<edge_reverse_t, graph_traits<BaseGraph>::edge_descriptor >> >> > Graph;

	size_t n = gt.size();
	size_t m = detected.size();
	size_t predStart = n + 1;
	size_t S = n, T = n + m + 1;

	std::vector<decltype(computeBand(detected[0], proximityTolerance))> bands;
	for (size_t i = 0; i < m; ++i)
		bands.push_back(computeBand(detected[i], proximityTolerance));

	Graph g(n + m + 2);

	// add S-edges
	for (int i = 0; i < n; ++i) {
		auto e = add_edge(S, i, g);
		auto er = add_edge(i, S, g);
		put(edge_reverse_t(), g, e.first, er.first);
		put(edge_reverse_t(), g, er.first, e.first);
		put(edge_capacity_t(), g, e.first, 1);
	}
	// add T-edges
	for (size_t i = 0; i < m; ++i) {
		auto e = add_edge(i + predStart, T, g);
		auto er = add_edge(T, i + predStart, g);
		put(edge_reverse_t(), g, e.first, er.first);
		put(edge_reverse_t(), g, er.first, e.first);
		put(edge_capacity_t(), g, e.first, 1);
	}
	// main edges
	int edgeCount = 0;
	for (size_t i = 0; i < S; ++i)
		for (size_t j = predStart; j < T; ++j) {
			float score = computeScore(gt[i], bands[j - predStart]);
			if (score) {
				auto e = add_edge(i, j, g);
				auto er = add_edge(j, i, g);
				put(edge_reverse_t(), g, e.first, er.first);
				put(edge_reverse_t(), g, er.first, e.first);
				put(edge_capacity_t(), g, e.first, 1);
				put(edge_weight_t(), g, e.first, 1 - score);
				put(edge_weight_t(), g, er.first, score - 1);
				++edgeCount;
			}
		}
	
	std::cout << edgeCount << " edges added, " << n << " hairs in GT, " << m << " tracks detected" << std::endl;

	// solve
	successive_shortest_path_nonnegative_weights(g, S, T);

	// read matching
	graph_traits<Graph>::edge_iterator ei, ei_end;
	std::vector<double> coverage(10);
	double total = 0, totalRelevant = 0;
	int matchedHairs = 0;
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		int res = get(edge_residual_capacity_t(), g, *ei);
		if (res == 0 && source(*ei, g) < S && predStart <= target(*ei, g) && target(*ei, g) < T) {
			++matchedHairs;
			float weight = 1 - get(edge_weight_t(), g, *ei);
			total += weight;
			++coverage[weight < 1 ? size_t(weight * coverage.size()) : coverage.size() - 1];
			totalRelevant += computeScore(detected[target(*ei, g) - predStart], computeBand(gt[source(*ei, g)], proximityTolerance));
		}
	}
	std::cout 
		<< matchedHairs << " hairs matched" << std::endl
		<< "total coverage (recall): " << total / n << " (" << total / matchedHairs << ")" << std::endl
		<< "precision: " << totalRelevant / m << std::endl;
	std::cout << "coverage histogram: " << std::endl;

	printf("       ");
	for (size_t i = 0; i < coverage.size(); ++i) {
		printf("%7d%%", (i + 1) * 100 / coverage.size());
	}
	printf("%8s", "Total");
	std::cout << std::endl;

	printf("cnt    "); double covSum = 0;
	for (size_t i = 0; i < coverage.size(); ++i) {
		printf("%8.0lf", coverage[i]);
		covSum += coverage[i];
	}
	printf("%8.0lf", covSum);
	std::cout << std::endl;

	printf("%%      "); covSum = 0;
	for (size_t i = 0; i < coverage.size(); ++i) {
		printf("%8.2lf", coverage[i] / n * 100);
		covSum += coverage[i] / n * 100;
	}
	printf("%8.2lf", covSum);
	std::cout << std::endl;

}

int main(int argc, const char* argv[]) {

	//evaluate({ { { 0.f, 0 }, { 1, 1 } } }, { { { 0.f, 1 }, { 1, 0 }, { 2, 1 } } }, 00.5);

	std::string inputPath;
	std::string calibPath;
	std::string imagePath;
	std::string resultFolder;
	std::string maskPath;

	std::string groundTruthPath;
	float proximityTolerance;

	cv::Size size;

	po::options_description desc("Allowed options");
	desc.add_options()
		("help", "produce help message")
		("i", po::value(&inputPath)->required(), "input *.txt_3dlines")
		("c", po::value(&calibPath)->required(), "calibration parameters file")
		("gt", po::value(&groundTruthPath)->required(), "ground truth path template *.txt_3dlines")
		("mask", po::value(&maskPath)->required(), "mask path template *.png")
		("o", po::value(&resultFolder)->required(), "destination folder")
		("tol", po::value(&proximityTolerance)->default_value(2), "tolerance in px")
		;

	po::positional_options_description p;
	p.add("i", 1);
	p.add("c", 1);
	p.add("gt", 1);
	p.add("o", 1);

	po::variables_map vm;

	bool testRun = true;

	try {
		po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
		po::notify(vm);

		if (vm.count("help")) {
			std::cout << desc << "\n";
			return 1;
		}

		CCameraArray cameras(calibPath.c_str());

		std::vector<std::vector<CVec3f>> tracks;
		std::ifstream ifs(inputPath);
		std::string tmp;
		for (std::vector<CVec3f> track; 
			readVector(ifs, track); 
			std::getline(ifs, tmp), tracks.push_back(track));

		std::vector<std::vector<CVec2f>> detected;

		for (unsigned c = 0; c < cameras.ArrayNumCams(); ++c) {
			CCamera& cam = cameras.GetNthCam(c);

			std::string gtPath = format(groundTruthPath.c_str(), cam.Id);
			if (!boost::filesystem::exists(gtPath))
				continue;

			for (int i = 0; i < tracks.size(); ++i) {
				detected.push_back(project(cam, tracks[i]));
			}


			std::vector<std::vector<CVec2f>> gt;
			std::ifstream ifs(gtPath);
			for (std::vector<CVec2f> track; readVector(ifs, track); gt.push_back(track));

			if (testRun) {
				std::cout << std::endl << std::endl << "Test run..." << std::endl;
				evaluate(gt, gt, proximityTolerance);
				testRun = false;
			}

			std::cout << std::endl << std::endl << "Evaluation againt ground truth for image #" << cam.Id << std::endl;

			if (vm.count("mask")) {
				std::cout << "Filtering GT tracks" << std::endl;

				std::string path = format(maskPath.c_str(), cam.Id);
				CV_Assert(boost::filesystem::exists(path));
				Mat mask = imread(path, CV_LOAD_IMAGE_UNCHANGED);

				auto bad = [&mask](const std::vector<CVec2f>& track) {
					for (size_t i = 1; i < track.size(); ++i) {
						LineIterator it(mask, vector_cast<Point>(track[i - 1]), vector_cast<Point>(track[i]));
						for (int i = 0; i < it.count; ++i, ++it)
							if (**it)
								return false;
					}
					return true;
				};
				gt.erase(std::remove_if(gt.begin(), gt.end(), bad), gt.end());
			}

			evaluate(detected, gt, proximityTolerance);
		}
	}
	catch (const std::exception& ex) {
		std::cerr << ex.what() << std::endl << desc << std::endl;
		return -1;
	}
}
