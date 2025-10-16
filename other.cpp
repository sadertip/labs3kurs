#pragma once
#include "Matrices_template.h"
#include"Gaus.h"
#include "other.h"
#include <string>

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
	T ma = 0;
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
T norm1(const matrix<T>& A)
{
	T ma = 0;
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
	T ma = 0;
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
T norminf(const matrix<T>& A)
{
	T ma = 0;
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
matrix<T> b_mod(matrix<T>& b)
{
	auto res = b;
	for (int i = 1; i < res.cols * res.rows; i += 2)
	{
		res[i - 1] += 0.01;
		res[i] -= 0.01;
	}
	return res;
}

template<typename T>
void matrix_test(const std::string& matrix_file, const std::string& vector_file)
{
	matrix<T> A = matrix<T>(matrix_file);
	matrix<T> b = matrix<T>(vector_file);
	std::cout << "Matrix A: \n" << A << '\n' << "Vector b: \n" << b << std::endl;
	auto x1 = Gauss(A, b);
	auto x2 = QR_decomposition(A, b);

	std::cout << "Gauss ans: \n" << x1 << '\n' << "QR ans: \n" << x2 << std::endl;
	auto b1 = A * x1;
	auto residual1 = mismatch(b1, b);
	std::cout << "Ans residual = " << residual1 << '\n' << std::endl;

	matrix<T> b0 = b_mod(b);
	std::cout << "b mod: \n" << b0 << std::endl;

	matrix<T> x01 = Gauss(A, b0);

	auto ans_diff1 = mismatch(x01, x1);
	std::cout << "Ans mod: \n" << x01 << "Ans diff = " << ans_diff1 << '\n' << std::endl;

	auto condV1 = cond1(A);
	auto condVinf = condinf(A);
	std::cout << "Matrix cond1 = " << condV1 << "\nMatrix condinf = " << condVinf << std::endl;


	//auto b0 = matrix<T>("..\\labs3kurs\\tests\\test1_vector1.dat");


	auto Dx = mismatch(x1, x01);
	auto dx2 = Dx / norm2(x1);
	//auto dxinf = Dx / norminf(x1);

	auto Db = mismatch(b, b0);
	auto db2 = Db / norm2(b0);
	//auto dbinf = Db / norminf(x1);

	auto cond_est2 = dx2 / db2;
	//auto cond_estinf = dxinf / dbinf;

	std::cout << "Matrix cond estimation = " << cond_est2 << std::endl;
	//std::cout << "\nMatrix cond estimation with norminf = " << cond_estinf << std::endl;

	auto A1 = reverse(A);
	auto E = A1 * A;
	std::cout << "\nE: \n" << E << std::endl;
}

template<typename T>
T mismatch(matrix<T>& a, matrix<T>& b)
{
	T sum = 0;
	for (size_t i = 0; i < a.rows; ++i)
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

		for (size_t i = 0; i < n; ++i) {
			b_vector(i, 0) = (i == j) ? T(1) : T(0);
		}

		matrix<T> A_copy = m;
		matrix<T> solution_vector = Gauss(A_copy, b_vector);

		for (size_t i = 0; i < n; ++i) {
			result_matrix(i, j) = solution_vector(i, 0);
		}
	}

	return result_matrix;
}

template<typename T>
T cond1(matrix<T>& m)
{
	matrix<T> m1 = reverse(m);
	T res = norm1(m1) * norm1(m);
	return res;
}

template<typename T>
T cond2(matrix<T>& m)
{
	matrix<T> m1 = reverse(m);
	T res = norm2(m1) * norm2(m);
	return res;
}

template<typename T>
T condinf(matrix<T>& m)
{
	matrix<T> m1 = reverse(m);
	T res = norminf(m1) * norminf(m);
	return res;
}