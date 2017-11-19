#pragma once

#include "types.h"
#include "hermite_spline.h"
#include "ceres.h"
#include "geom.h"

#include <Eigen/Dense>

#include <ceres/ceres.h>
#include <ceres/rotation.h>

struct HermiteEnergy {
	const Track& target;
	double maxLen;
	double damping;

	template<class T>
	bool operator() (const T* params, T* output) const {
		typedef Eigen::Matrix<T, 3, 1> Vec;
		typedef Eigen::Map<Vec> Map;
		typedef Eigen::Map<const Vec> CMap;

		Vec tBegin = target.front().template cast<T>();
		Vec tEnd = target.back().template cast<T>();

		auto delta1 = (target[1] - target[0]).template cast<T>();
		CMap m1(params);
		Map(output + 0) = T(damping) *  m1.cross(delta1) / sqrt(m1.dot(m1) * delta1.dot(delta1));
		auto delta2 = (target[target.size() - 2] - target[target.size() - 1]).template cast<T>();
		CMap m2(params + 3);
		Map(output + 3) = T(damping) * m2.cross(delta2) / sqrt(m2.dot(m2) * delta2.dot(delta2));

		std::vector<Vec> track, ref;
		for (int i = 0; i < 100; ++i)
			track.push_back(hermite_spline::val(CMap(tBegin.data()), m1, CMap(tEnd.data()), m2,
				hermite_spline::H(i / 99.)));
		for (auto& p : target)
			ref.emplace_back(p.cast<T>());
		for (size_t j = 1; j + 1 < target.size(); ++j) {
			double bestDistance = INFINITY;
			Vec diff;
			distance(ref[j].data(), track, diff.data());
			Map(output + 3 * j + 3) = diff;
		}
		return true;
	}
};

ceres::Solver::Summary fitSplineToTrack(const Track& track, double maxLen, double damping, double parameters[6]) {

	Map m0(parameters + 0);
	m0 = track[1] - track[0];
	if (maxLen / m0.norm() > .5)
		m0 *= maxLen * 2 / m0.norm();

	Map m1(parameters + 3);
	m1 = track[track.size() - 1] - track[track.size() - 2];
	if (maxLen / m1.norm() > .5)
		m1 *= maxLen * 2 / m1.norm();

	ceres::Problem problem;
	problem.AddResidualBlock(new ceres::AutoDiffCostFunction<HermiteEnergy, ceres::DYNAMIC, 6>(
		new HermiteEnergy{ track, maxLen, damping }, track.size() * 3),
		NULL,
		parameters);

	ceres::Solver::Options options;
	options.max_num_iterations = 1000;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	return summary;
}
