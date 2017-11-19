#pragma once

#include <Eigen/Core>

#include <ios>

template<class T, int R, int C>
std::istream& operator>>(std::istream& is, Eigen::Matrix<T, R, C>& mat) {
	for (int k = 0; k < R*C; ++k) {
		is >> mat.data()[k];
	}
	return is;
}

template<class T, int R, int C>
std::ostream& operator<<(std::ostream& os, Eigen::Matrix<T, R, C>& mat) {
	for (int k = 0; k < R*C; ++k) {
		os << mat.data()[k] << " ";
	}
	return os;
}
