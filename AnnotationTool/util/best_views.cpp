#include "best_views.h"

#include "../eval_reconstruction/successive_shortest_path_nonnegative_weights.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <opencv2/core.hpp>

#include <iostream>
#include <vector>

std::vector<std::vector<unsigned>> determineBestViews(const std::vector<std::vector<double>>& scores, size_t topK, double viewChangeCost, size_t* changes, double* totalCost) {
	double absMin = INFINITY;
	for (auto& vec : scores) {
		CV_Assert(vec.size() == scores[0].size());
		for (auto& val : vec)
			absMin = std::min(absMin, val);
	}

	using namespace boost;

	typedef adjacency_list<vecS, vecS, directedS> BaseGraph;

	typedef adjacency_list<vecS, vecS, directedS, no_property,
		property<edge_weight_t, double,
		property<edge_capacity_t, int,
		property<edge_residual_capacity_t, int,
		property<edge_reverse_t, graph_traits<BaseGraph>::edge_descriptor >> >> > Graph;

	size_t viewCount = scores[0].size();
	size_t segmentCount = scores.size();
	size_t nodes = viewCount * segmentCount * 2 + 3;
	size_t dualShift = viewCount * segmentCount;

	Graph g(nodes);

	size_t S = nodes - 3,
		S1 = nodes - 2,
		T = nodes - 1;

	// add S-edges
	auto e = add_edge(S, S1, g);
	auto er = add_edge(S1, S, g);
	put(edge_reverse, g, e.first, er.first);
	put(edge_reverse, g, er.first, e.first);
	put(edge_capacity, g, e.first, topK);

	// add S1-edges
	for (int i = 0; i < viewCount; ++i) {
		e = add_edge(S1, i, g);
		er = add_edge(i, S1, g);
		put(edge_reverse, g, e.first, er.first);
		put(edge_reverse, g, er.first, e.first);
		put(edge_capacity, g, e.first, 1);
	}

	// add T-edges
	for (size_t i = 0, shift = (segmentCount - 1) * viewCount + dualShift; i < viewCount; ++i) {
		e = add_edge(i + shift, T, g);
		er = add_edge(T, i + shift, g);
		put(edge_reverse, g, e.first, er.first);
		put(edge_reverse, g, er.first, e.first);
		put(edge_capacity, g, e.first, 1);
	}

	// add inter-edges
	for (size_t i = 0; i < segmentCount; ++i)
		for (size_t j = 0; j < viewCount; ++j) {
			size_t id = i * viewCount + j;
			e = add_edge(id, id + dualShift, g);
			er = add_edge(id + dualShift, id, g);
			put(edge_reverse, g, e.first, er.first);
			put(edge_reverse, g, er.first, e.first);
			put(edge_capacity, g, e.first, 1);
			put(edge_weight, g, e.first, scores[i][j] - absMin);
			put(edge_weight, g, er.first, absMin - scores[i][j]);
		}

	// add view-transition edges
	for (size_t i = 0; i + 1 < segmentCount; ++i)
		for (size_t k = 0; k < viewCount; ++k) {
			for (size_t j = 0; j < viewCount; ++j) {
				size_t idj = (i + 1) * viewCount + j;
				size_t idk = i * viewCount + k;

				e = add_edge(idk + dualShift, idj, g);
				er = add_edge(idj, idk + dualShift, g);
				put(edge_reverse, g, e.first, er.first);
				put(edge_reverse, g, er.first, e.first);
				put(edge_capacity, g, e.first, 1);
				if (k != j) {
					put(edge_weight, g, e.first, viewChangeCost);
					put(edge_weight, g, er.first, -viewChangeCost);
				}
			}
		}

	// solve
	successive_shortest_path_nonnegative_weights(g, S, T);

	// read matching
	graph_traits<Graph>::edge_iterator ei, ei_end;
	std::vector<std::vector<unsigned>> result(segmentCount);
	size_t viewChanges = 0;
	for (boost::tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
		int res = get(edge_residual_capacity, g, *ei);
		int cap = get(edge_capacity, g, *ei);
		size_t src = source(*ei, g);
		size_t trg = target(*ei, g);
		if (res == 0 && cap > 0) { // if there is flow 
			if (src < dualShift && src + dualShift == trg) {
				size_t segment = src / viewCount;
				size_t view = src % viewCount;
				result[segment].push_back(view);
			}
			if (dualShift <= src && src < S && trg < dualShift) {
				viewChanges += src - dualShift != trg - viewCount;
			}
		}
	}

#ifdef _DEBUG
	std::cout << "Total changes: " << viewChanges << std::endl;
	std::cout << "Segment-to-segment average changes: " << (double)viewChanges / segmentCount << std::endl;
#endif

	if (changes)
		*changes = viewChanges;
	if (totalCost)
		*totalCost = boost::find_flow_cost(g, /* avoid bug https://svn.boost.org/trac/boost/ticket/12700 */ bgl_named_params<int, buffer_param_t>())
			+ absMin * topK * segmentCount;

	// done!
	return result;
}