#pragma once
#include "Matrices_template.h"
#include"Gaus.h"
#include "other.h"
template<typename T>
T norm2(matrix<T>& A)
{
	T sum_of_squares = 0;
	for (size_t i = 0; i < A.rows; ++i) {
		T sumi = 0;
		for (size_t j = 0; j < A.cols; ++j) {
			sumi += A(i, j) * A(i, j);
		}
		sum_of_squares += sumi;
	}
	return std::sqrt(sum_of_squares);
}
template<typename T>
T norm1(matrix<T>& A)
{
	T ma = A(0, 0);
	for (size_t j = 0; j < A.cols; ++j) {
		T sumi = 0;
		for (size_t i = 0; i < A.rows; ++i) {
			sumi += std::abs(A(i, j));
		}
		ma = std::max(ma, sumi);
	}
	return ma;
}
template<typename T>
T norminf(matrix<T>& A)
{
	T ma = A(0, 0);
	for (size_t i = 0; i < A.rows; ++i) {
		T sumi = 0;
		for (size_t j = 0; j < A.cols; ++j) {
			sumi += std::abs(A(i, j));
		}
		ma = std::max(ma, sumi);
	}
	return ma;
}

template<typename T>
T mismatch(matrix<T>& a, matrix<T>& b)
{
	T sum = 0;
	for (size_t i = 0; i <= a.rows; ++i)
	{
		sum += (a[i] - b[i]) * (a[i] - b[i]);
	}
	return std::sqrt(sum);
}

template<typename T>
matrix<T> reverse(const matrix<T>& m) {
	if (m.rows != m.cols) {
		throw std::invalid_argument("Matrix must be square to be invertible.");
	}
	size_t n = m.rows;

	matrix<T> result_matrix(n, n);

	matrix<T> b_vector(n, 1);

	for (size_t j = 0; j < n; ++j) {
		// Формирование вектора b (e_j)
		for (size_t i = 0; i < n; ++i) {
			b_vector(i, 0) = (i == j) ? T(1) : T(0);
		}

		matrix<T> A_copy = m;


		matrix<T> solution_vector = Gaus(A_copy, b_vector);

		for (size_t i = 0; i < n; ++i) {
			result_matrix(i, j) = solution_vector(i, 0);
		}
	}

	return result_matrix;
}

template<typename T>
T cond(matrix<T>& m)
{
	matrix<T> m1 = reverse(m);
	T res = norm2(m1) * norm2(m);
	return res;
}