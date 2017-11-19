#pragma once

#include "ceres_fwd.h"

#include <ceres/ceres.h>

#include <opencv2/core.hpp>

template<class T, int N>
struct GetType<ceres::Jet<T, N>> {
	typedef T type;
	static T removeJet(const ceres::Jet<T, N>& t) { return t.a; }
};

template<class T>
void ok(const T& t) { CV_Assert(std::isnormal(t) || t == 0); }

template<class T, int N>
void ok(const ceres::Jet<T, N>& t) {
	ok(t.a);
	for (int i = 0; i < N; ++i)
		ok(t.v[i]);
}

template<class T, int N>
void ok(const Eigen::Matrix<T, N, 1>& t) {
	for (int i = 0; i < N; ++i)
		ok(t[i]);
}

template<class T, int N>
void ok(const T (&t)[N]) {
	for (int i = 0; i < N; ++i)
		ok(t[i]);
}

template<class R, class T, int N>
R saturate_cast(const ceres::Jet<T, N>& jet) {
	using cv::saturate_cast;
	return saturate_cast<R>(jet.a);
}

template<int N, class T>
Eigen::Map<const Eigen::Matrix<T, N, 1>> removeJet(const T* params) {
	return Eigen::Map<const Eigen::Matrix<T, N, 1>> (params);
}

template<int N, class T, int M>
Eigen::Matrix<T, N, 1> removeJet(const ceres::Jet<T, M>* params) {
	Eigen::Matrix<T, N, 1> res;
	for (int i = 0; i < N; ++i)
		res[i] = asConst(params[i]);
	return res;
}

template<class Functor, int Params>
class GradientFunction : public ceres::FirstOrderFunction {
public:
	Functor* func;

	explicit GradientFunction(Functor* f = new Functor())
		: func(f)
		, autodiff(f)
	{
	}

	bool Evaluate(const double*const parameters,
		double* cost,
		double* gradient) const override
	{
		return autodiff.Evaluate(&parameters, cost, &gradient);
	}

	int NumParameters() const override { return Params; }

private:
	ceres::AutoDiffCostFunction<Functor, 1, Params> autodiff;
};

template<class T, class Mat34>
void project(const Mat34& mat, const T p[3], T v[2], T* depth = NULL) {
	typedef typename GetType<T>::type Real;
	v[1] = Real(mat(2, 0)) * p[0] + Real(mat(2, 1)) * p[1] + Real(mat(2, 2)) * p[2] + Real(mat(2, 3));
	if (depth)
		*depth = v[1];
	v[0] = (Real(mat(0, 0)) * p[0] + Real(mat(0, 1)) * p[1] + Real(mat(0, 2)) * p[2] + Real(mat(0, 3))) / v[1];
	v[1] = (Real(mat(1, 0)) * p[0] + Real(mat(1, 1)) * p[1] + Real(mat(1, 2)) * p[2] + Real(mat(1, 3))) / v[1];
}

template<class DataType, class Derived>
typename Derived::Scalar bilinear(const cv::Mat_<DataType>& data, const Eigen::MatrixBase<Derived>& p) {
	typedef typename Derived::Scalar T;

	cv::Rect box(cv::Point(), data.size());
	cv::Point currentPos((int)asConst(p[0]), (int)asConst(p[1]));

	//bilinear interpolation
	T dx = p[0] - double(currentPos.x);
	T dy = p[1] - double(currentPos.y);

	auto getVal = [&](int y, int x) {
		cv::Point p = currentPos + cv::Point(x, y);
		if (box.contains(p)) {
			return double(data(p));
		}
		return double(0);
	};

	T value =
		(1.0 - dx) * (1.0 - dy) * getVal(0, 0)
		+ (dx)* (1.0 - dy)       * getVal(0, 1)
		+ (1.0 - dx) * (dy)* getVal(1, 0)
		+ (dx)* (dy)* getVal(1, 1);

	return value;
}
