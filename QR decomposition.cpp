#pragma once
#include <iostream>
#include "QR decomposition.h"
//#include "Matrices_template.cpp"
#include "Gaus.cpp"

template<typename T>
matrix<T> QR_decomposition(matrix<T>& A, matrix<T>& b)
{
	matrix<T> Q = eye<T>(A.rows, 1);
	matrix<T> R = matrix(A);
	matrix<T> b1 = matrix(b);

	for (int step = 0; step < A.cols - 1; ++step)
	{
		QR_col_step(R, b1, Q, step);
	}

	auto x = gauss_back(R, b1);
	std::cout << "Matrix Q = \n" << Q << "Matrix R = \n" << R << std::endl;
	//auto A1 = Q * R;
	//std::cout << A1 << std::endl;
	return x;
}

template<typename T>
std::pair<T, T> Coeffs_ij(matrix<T>& A, int i, int j)
{
	T a1 = A(i, i);
	T a2 = A(j, i);
	T radix = sqrt(a1 * a1 + a2 * a2);
	T c = a1 / radix;
	T s = a2 / radix;
	return { c,s };
}

template<typename T>
void QR_row_step(matrix<T>& A, matrix<T>& b, T c_ij, T s_ij, int row, int step)
{
	T* temp = new T[A.cols];

	//std::cout << "before: \n" << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	for (int j = 0; j < A.cols; ++j)
	{
		temp[j] = A[A.pos(step, j)];
	}
	T b_temp = b[step];

	for (int j = 0; j < A.cols; ++j)
	{
		A[A.pos(step, j)] = c_ij * A[A.pos(step, j)] + s_ij * A[A.pos(row, j)];
	}
	b[step] = c_ij * b_temp + s_ij * b[row];

	//std::cout << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	for (int j = 0; j < A.cols; ++j)
	{
		A[A.pos(row, j)] = -s_ij * temp[j] + c_ij * A[A.pos(row, j)];
	}
	b[row] = -s_ij * b_temp + c_ij * b[row];

	//std::cout << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	delete[] temp;
}

template<typename T>
void QR_col_step(matrix<T>& A, matrix<T>& b, matrix<T>& Q, int step)
{
	std::pair<T, T> coeffs;

	for (int i = step + 1; i < A.rows; ++i)
	{
		coeffs = Coeffs_ij(A, step, i);

		matrix<T> Temp = eye<T>(Q.rows, 1.);

		Temp(step, step) = coeffs.first;
		Temp(step, i) = coeffs.second;
		Temp(i, i) = coeffs.first;
		Temp(i, step) = -coeffs.second;

		//std::cout << "step, row = " << step << ' ' << i << '\n' << Temp << std::endl;

		Q = Temp * Q;

		//std::cout << Q << std::endl;

		QR_row_step(A, b, coeffs.first, coeffs.second, i, step);
	}
}

//int main()
//{
//	using T = double;
//	matrix<T> A = matrix<T>(3, 3, { 1,0,4,0,1,1,1,2,1 });
//	matrix<T> b = matrix<T>(3, 1, { 1, 2, 3 });
//	QR_decomposition(A, b);
//	return 0;
//}