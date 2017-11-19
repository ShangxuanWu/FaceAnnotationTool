#pragma once

#include <utility>

#include <Eigen/Dense>

template<class T>
struct GetType {
	typedef T type;

	static type removeJet(const T& t) {
		return t;
	}
};

template<class T, int N, int M>
struct GetType<Eigen::Matrix<T, N, M>> {
	typedef Eigen::Matrix<typename GetType<T>::type, N, M> type;
	static type removeJet(const Eigen::Matrix<T, N, M>& t)  {
		return t.unaryExpr(std::ptr_fun(GetType<T>::removeJet));
	}
};

template<class T>
typename GetType<T>::type asConst(const T& t) { return GetType<T>::removeJet(t); }
