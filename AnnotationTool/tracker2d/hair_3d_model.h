#pragma once

#include "../detect_1view/pca_model.h"
#include "../detect_1view/rotation_multi.h"

#include "../util/ceres.h"
#include "../util/geom.h"

#include <random>

template<int N>
struct Hair3DModel {

	typedef RotationMulti_Model<PCANormal_Model<N>> Model;

	struct ProjectionEnergy {
		const Model& model;
		const TrackT<double, 2>& target;
		const CMtx4x4f& P;
		double w;

		template<class T>
		bool operator()(const T* params, T* output) const {
			auto track = generateHair(model, params);
			TrackT<T, 2> projection, targetT;
			for (auto& p : track) {
				VecT<T, 2> q;
				project(P, p.data(), q.data());
				projection.push_back(q);
			}
			for (auto& p : target)
				targetT.push_back(p.template cast<T>());
			
			//computeDistances(projection, targetT, output + targetT.size() * 3);
			computeDistances(targetT, projection, output);
			return true;
		}

		template<class T>
		void computeDistances(const TrackT<T, 2>& from, const TrackT<T, 2>& to, T* output) const {
			T f = T(w / sqrt((double)from.size()) / 2);
			for (size_t i = 0; i < from.size(); ++i) {
				auto& p = from[i];
				MapT<T, 2> q(output + i * 3);
				distance<2>(p.data(), to, q.data());
				q *= f;

				if (from.size() > 1) {
					VecT<T, 2> norm;
					if (i == 0)
						norm = from[1] - from.front();
					else if (i == from.size() - 1)
						norm = from[i] - from[i - 1];
					else
						norm = from[i + 1] - from[i - 1];
					std::swap(norm(0), norm(1));
					norm(1) = -norm(1);
					output[i * 3 + 2] = q.dot(norm.normalized()) * f;
				}
				else {
					output[i * 3 + 2] = T(0);
				}
			}
		}
	};

	struct ModelEnergy {
		const Model& model;

		template<class T>
		bool operator()(const T* params, T* output) const {
			return computeCost(model, params, output);
		}
	};

	Model model;
	double params[Model::VAR_PARAMS_NUMBER];
	double curParams[Model::VAR_PARAMS_NUMBER];
	CMtx4x4f P;
	double projWeight;
	std::mt19937 gen;

	Hair3DModel(const std::string& clusters, int target, double w, 
		const std::string& meanPath, const std::string& eigenPath, double pcaW, 
		const std::string& scalePath, double scaleW,
		const CMtx4x4f& P, double projWeight)
		: model(clusters, target, w, meanPath, eigenPath, pcaW, scalePath, scaleW)
		, P(P)
		, projWeight(projWeight)
	{
		reset();
	}

	void reset() {
		sampleHair(model, params, gen);
	}

	double adjust(const TrackT<double, 2>& projection, double curParams[]) const {
		memcpy(curParams, params, sizeof(params));

		ceres::Problem problem;
		problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ModelEnergy, Model::OUTPUT_NUMBER, Model::VAR_PARAMS_NUMBER>(
			new ModelEnergy{ model }), NULL, curParams);
		
		problem.AddResidualBlock(new ceres::AutoDiffCostFunction<ProjectionEnergy, ceres::DYNAMIC, Model::VAR_PARAMS_NUMBER>(
			new ProjectionEnergy{ model, projection, P, projWeight }, projection.size() * 3), NULL, curParams);

		ceres::Solver::Summary summary;
		ceres::Solver::Options options;
		//options.minimizer_progress_to_stdout = true;
		options.max_num_iterations = 1000;

		ceres::Solve(options, &problem, &summary);


		//std::cout << summary.BriefReport() << std::endl; 

		return summary.final_cost;
	}

	void commit(double curParams[]) {
		memcpy(params, curParams, sizeof(params));
	}
};
