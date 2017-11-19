#pragma once

#include "util.h"
#include "ceres.h"

#include <Eigen/Core>

#include <algorithm>

namespace hermite_spline {
	template<class T, int N = 3>
	using Vec = Eigen::Matrix<T, N, 1>;

	template<class T>
	static Vec<T, 4> H(const T& t) {
		return{ (1. + 2. * t)*sqr(1. - t), t*sqr(1. - t), sqr(t)*(3. - 2. * t), sqr(t)*(t - 1.) };
	}

	template<class T>
	static Vec<T, 4> Hprime(const T& t) {
		return{ -6. * t*(1. - t), (1. - t)*(1. - 3. * t), 6. * t*(1. - t), t*(3. * t - 2.) };
	}

	template<class T>
	static Vec<T, 4> Hpp(const T& t) {
		return{ 12. * t - 6., 6. * t - 4., 6. - 12. * t, 6. * t - 2. };
	}

	template<class Derived, class R>
	static auto val(
		const Eigen::MatrixBase<Derived>& begin,
		const Eigen::MatrixBase<Derived>& m0,
		const Eigen::MatrixBase<Derived>& end,
		const Eigen::MatrixBase<Derived>& m1,
		const Vec<R, 4>& base)

		-> decltype(begin(0) * begin + begin(0) * begin + begin(0) * begin + begin(0) * begin)
	{
		typedef typename Derived::Scalar T;
		return T(base[0]) * begin + T(base[2]) * end + T(base[1]) * m0 + T(base[3]) * m1;
	}

	template<class T, class Visitor>
	static Visitor& visit(const T* _begin, const T* _m0, const T* _end, const T* _m1, double maxLen, Visitor& visitor) {
#define AS_EIGEN_VEC(name) Vec<T> name(_ ## name); Vec<double> name ## const(asConst(name(0)), asConst(name(1)), asConst(name(2)))
		AS_EIGEN_VEC(begin);
		AS_EIGEN_VEC(m0);
		AS_EIGEN_VEC(end);
		AS_EIGEN_VEC(m1);
		for (double t = 0;;) {
			auto cur = val(begin, m0, end, m1, H(t));

			visitor(begin, m0, end, m1, t, cur);

			if (t >= 1)
				break;
			double step;
			auto length = val(beginconst, m0const, endconst, m1const, Hprime(t)).norm();
			if (length > 0)
				step = maxLen / length;
			else {
				length = val(beginconst, m0const, endconst, m1const, Hpp(t)).norm();
				if (length > 0)
					step = sqrt(2 * maxLen / length);
				else if (begin == end) {
					std::cerr << "begin == end" << std::endl;
					return visitor;
				}
				else {
					step = 1e-6;
					std::cerr << "Two derivatives are zero!!!" << std::endl;
				}
			}
			step = std::min(step, 0.5);
			CV_Assert(step > 0);
			if (cur == val(begin, m0, end, m1, H(std::min(1.0, t + step))))
				std::cout << 1 - t << " " << t + step << std::endl;
			t = std::min(1.0, t + step);
			if (1 - t < step)
				t = 1;
		}

		return visitor;
	}

	template<class T, class It>
	It rasterize(const T* _begin, const T* _m0, const T* _end, const T* _m1, double maxLen, It out) {
		struct Visitor {
			void operator()(const Vec<T>&, const Vec<T>&, const Vec<T>&, const Vec<T>&, double, const Vec<T>& cur) {
				*out++ = cur;
			}
			It out;
		} visitor{ out };
		return visit(_begin, _m0, _end, _m1, maxLen, visitor).out;
	}
};
