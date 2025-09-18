#pragma once
#include <iostream>
#include "QR decomposition.h"
//#include "Matrices_template.cpp"
#include "Gaus.cpp"

matrix<double> QR_decomposition(matrix<double>& A, matrix<double>& b)
{
	matrix<double> Q = eye<double>(A.rows, 1);
	matrix<double> R = matrix<double>(A);
	matrix<double> b1 = matrix<double>(b);

	for (int step = 0; step < A.cols - 1; ++step)
	{
		QR_col_step(R, b1, Q, step);
	}

	auto x = gauss_back(R, b1);
	std::cout << "Matrix Q = \n" << Q << "Matrix R = \n" << R << "Ans = " << x << std::endl;
	auto A1 = Q * R;
	std::cout << A1 << std::endl;
	return x;
}

std::pair<double, double> Coeffs_ij(matrix<double>& A, int i, int j)
{
	double a1 = A[A.pos(i, i)];
	double a2 = A[A.pos(j, i)];
	double radix = sqrt(a1 * a1 + a2 * a2);
	double c = a1 / radix;
	double s = a2 / radix;
	return { c,s };
}

void QR_row_step(matrix<double>& A, matrix<double>& b, double c_ij, double s_ij, int row, int step)
{
	double* temp = new double[A.cols];

	std::cout << "before: \n" << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	for (int j = 0; j < A.cols; ++j)
	{
		temp[j] = A[A.pos(step, j)];
	}
	double b_temp = b[step];

	for (int j = 0; j < A.cols; ++j)
	{
		A[A.pos(step, j)] = c_ij * A[A.pos(step, j)] + s_ij * A[A.pos(row, j)];
	}
	b[step] = c_ij * b_temp + s_ij * b[row];

	std::cout << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	for (int j = 0; j < A.cols; ++j)
	{
		A[A.pos(row, j)] = -s_ij * temp[j] + c_ij * A[A.pos(row, j)];
	}
	b[row] = -s_ij * b_temp + c_ij * b[row];

	std::cout << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	delete[] temp;
}

void QR_col_step(matrix<double>& A, matrix<double>& b, matrix<double>& Q, int step)
{
	std::pair<double, double> coeffs;

	for (int i = step + 1; i < A.rows; ++i)
	{
		coeffs = Coeffs_ij(A, step, i);

		matrix<double> T = eye<double>(Q.rows, 1);

		T[T.pos(step, step)] = coeffs.first;
		T[T.pos(step, i)] = coeffs.second;
		T[T.pos(i, i)] = coeffs.first;
		T[T.pos(i, step)] = -coeffs.second;

		std::cout << "step, row = " << step << ' ' << i << '\n' << T << std::endl;

		Q = T * Q;

		std::cout << Q << std::endl;

		QR_row_step(A, b, coeffs.first, coeffs.second, i, step);
	}
}

int main()
{
	using T = double;
	matrix<T> A = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
	matrix<T> b = matrix<T>(3, 1, { 1, 2, 3 });
	QR_decomposition(A, b);
	return 0;
}