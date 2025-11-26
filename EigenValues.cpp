#include <iostream>
#include "EigenValues.h"
#include "Matrices_template.cpp"
#include "QR decomposition.cpp"
#include "simple_iteration_method.cpp"
#include "other.cpp"

int operations = 0;

template<typename T>
void Hessenberg_form(matrix<T>& A)
{
	for (int i = 1; i < A.cols; ++i)
	{
		for (int j = i + 1; j < A.rows; ++j)
		{
			rotation(A, i, j);
		}
	}
}

template<typename T>
std::pair<T, T> Alpha_beta(matrix<T>& A, int k, int l)
{
	T a1 = A(k, k - 1);
	T a2 = A(l, k - 1);
	T radix = sqrt(a1 * a1 + a2 * a2);
	T alpha = a1 / radix;
	T beta = a2 / radix;
	operations += 5;
	return std::pair<T, T>(alpha, beta);
}

template<typename T>
void rotation(matrix<T>& A, int k, int l)
{
	matrix<T> R = eye<T>(A.cols, 1);
	std::pair<T, T> AlphaBeta = Alpha_beta(A, k, l);
	R(k, k) = AlphaBeta.first;
	R(l, l) = AlphaBeta.first;
	R(k, l) = -AlphaBeta.second;
	R(l, k) = AlphaBeta.second;
	A = A * R;
	R(k, l) = AlphaBeta.second;
	R(l, k) = -AlphaBeta.second;
	A = R * A;
	operations += 8;
}

template<typename T>
T sigma_value(matrix<T>& A)
{
	T sigma = A(A.cols - 1, A.cols - 1);
	return sigma;
}

template<typename T>
void shift(matrix<T>& A, T sigma)
{
	auto E = eye<T>(A.cols, sigma);
	A = A - E;
}

template<typename T>
void backshift(matrix<T>& A, T sigma)
{
	auto E = eye<T>(A.cols, sigma);
	A = A + E;
}

template<typename T>
std::pair<matrix<T>, matrix<T>> QR_matrices(matrix<T>& A, bool Hessenberg)
{
	matrix<T> Q = eye<T>(A.rows, 1);
	matrix<T> R = matrix<T>(A);
	for (int step = 0; step < A.cols - 1; ++step)
	{
		QR_matrices_col_step(R, Q, step, Hessenberg);
	}
	Q = Q.tranpose();
	return std::pair<matrix<T>, matrix<T>>(Q, R);
}

template<typename T>
void QR_matrices_row_step(matrix<T>& A, T c_ij, T s_ij, int row, int step, bool Hessenberg)
{
	T* temp = new T[A.cols];

	//std::cout << "before: \n" << "R = \n" << A << '\n' << "b1 = \n" << b << std::endl;

	for (int j = 0; j < A.cols; ++j)
	{
		temp[j] = A(step, j);
	}

	for (int j = step; j < A.cols; ++j)
	{
		A(step, j) = c_ij * A(step, j) + s_ij * A(row, j);
		operations += 2;
	}

	//std::cout << "R = \n" << A << '\n' << std::endl;

	for (int j = step; j < A.cols; ++j)
	{
		A(row, j) = -s_ij * temp[j] + c_ij * A(row, j);
		operations += 2;
	}

	//std::cout << "R = \n" << A << '\n' << std::endl;

	delete[] temp;
}

template<typename T>
void QR_matrices_col_step(matrix<T>& A, matrix<T>& Q, int step, bool Hessenberg)
{
	std::pair<T, T> coeffs;

	int iter = Hessenberg ? step + 2 : A.rows;
	for (int i = step + 1; i < iter; ++i)
	{
		coeffs = Coeffs_ij(A, step, i);
		operations += 5;

		matrix<T> Temp = eye<T>(Q.rows, 1.);

		Temp(step, step) = coeffs.first;
		Temp(step, i) = coeffs.second;
		Temp(i, i) = coeffs.first;
		Temp(i, step) = -coeffs.second;

		//std::cout << "step, row = " << step << ' ' << i << '\n' << Temp << std::endl;

		Q = Temp * Q;
		operations += 4 * Q.cols;
		//std::cout << Q << std::endl;

		QR_matrices_row_step(A, coeffs.first, coeffs.second, i, step, Hessenberg);
	}
}

template<typename T>
matrix<T> QR_eigenvalues(matrix<T>& Aorig, T epsilon, bool Hessenberg, bool shifting)
{
	matrix<T> A = Aorig;
	matrix<T> eigenvalues = matrix<T>(1, A.cols);
	int iterations = 0;
	if (Hessenberg)
	{
		Hessenberg_form(A);
	}
	//std::cout << "A: \n" << A << std::endl;
	for (int n = 0; n < eigenvalues.cols - 1; ++n)
	{
		for (int k = 0; k < 100; ++k)
		{
			iterations += 1;
			auto sigma = sigma_value(A);
			//sigma = 4;
			if (shifting)
			{
				shift(A, sigma);
			}
			//std::cout << "A: \n" << A << std::endl;
			auto QR = QR_matrices(A, Hessenberg);
			//std::cout << "Q: \n" << QR.first << std::endl;
			//std::cout << "R: \n" << QR.second << std::endl;
			//A = QR.first * QR.second;
			//std::cout << "A: \n" << A << std::endl;
			A = QR.second * QR.first;
			operations += A.rows * A.cols;
			//std::cout << "A: \n" << A << std::endl;
			if (shifting)
			{
				backshift(A, sigma);
			}
			//std::cout << "A: \n" << A << std::endl;
			auto flag = true;
			//std::cout << "A: \n" << A << std::endl;
			for (int t = 0; t < A.cols - 1; ++t)
			{
				if (abs(A(A.rows - 1, t)) > epsilon)
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				eigenvalues[n] = A(A.cols - 1, A.cols - 1);
				A = A.minor(A.cols - 1);
				break;
			}
		}
	}
	eigenvalues[eigenvalues.cols - 1] = A[0];
	std::cout << "Hessenberg = " << Hessenberg <<
		", shifting = " << shifting << ", iteration = " <<
		iterations << ", operations = " << operations << std::endl;
	operations = 0;
	return eigenvalues;
}

template<typename T>
matrix<T> reverse_iter_eigenvectors(matrix<T>& A, matrix<T>& eigenvalues, T epsilon)
{
	matrix<T> eigenvectors = matrix<T>(A.rows, A.cols);
	for (int i = 0; i < eigenvalues.cols; ++i)
	{
		matrix<T> x = matrix<T>(A.rows, 1, 0);
		x(0, 0) = 1;
		std::cout << x << std::endl;
		matrix<T> temp = eye(A.rows, -eigenvalues(0, i)) + A;
		std::cout << A << std::endl;
		std::cout << temp << std::endl;
		for (int j = 0; j < 100; ++j)
		{
			std::cout << x << std::endl;
			auto check = A * x - x * eigenvalues(0, i);

			std::cout << check << std::endl;
			if (norm2(check) < epsilon)
			{
				for (int k = 0; k < A.cols; ++k)
				{
					eigenvectors(k, i) = x(k, 0);
				}
				auto b = matrix<T>(A.rows, 1);

				break;
			}
			//const matrix<T> x0 = x;
			x = QR_decomposition(temp, x);
			x = x * (1 / norm2(x));
			//x = simple_iteration_method(temp, x0, x0, 0.001, 100, epsilon);
			//x = QR_decomposition(temp, x);
		}
	}
	return eigenvectors;
}

template<typename T>
std::pair<matrix<T>, matrix<T>> Rayleigh(matrix<T>& Aorig, T epsilon)
{
	auto A = Aorig;
	matrix<T> eigenvalues = matrix<T>(1, A.cols);
	matrix<T> eigenvectors = matrix<T>(A.rows, A.cols);
	for (int i = 0; i < eigenvalues.cols; ++i)
	{
		matrix<T> x = matrix<T>(A.cols, 1, 0);
		x(i, 0) = 1;
		for (int k = 0; k < 100; ++k)
		{
			std::cout << "x: \n" << x << std::endl;
			T rho = rho_value(A, x);
			std::cout << rho << std::endl;
			auto temp = eye(A.cols, -rho) + A;
			std::cout << temp << std::endl;
			auto check = A * x - x * rho;
			std::cout << check << std::endl;
			if (norm2(check) < epsilon)
			{
				eigenvalues[i] = (rho-1);
				for (int k = 0; k < A.cols; ++k)
				{
					eigenvectors(k, i) = x(k, 0);
				}
				//A = eye(A.cols, rho) - A;
				//auto b = b_vector(A);
				//auto b = matrix<T>(1, A.cols, 0);
				//b[b.cols - 1] = 1;
				//b[0] = (rho - x[x.rows - 1]) / x[0];
				//std::cout << b << std::endl;
				////x = x * (1 / (x[x.rows - 1]));
				//std::cout << "x: \n" << x << std::endl;
				//A = A - x * b;
				//A = reverse(A + eye<T>(A.cols, 1));
				//std::cout << "A: \n" << A << std::endl;
				//break;
			}
			x = QR_decomposition(A, x);
			x = x * (1 / norm2(x));
		}
		//for (int i = 1; i < eigenvalues.cols; ++i)
		//	eigenvalues[i] = eigenvalues
		//eigenvalues[0] = (1/eigenvalues[0]);
	}
	return std::pair<matrix<T>, matrix<T>>(eigenvalues,eigenvectors);
}

template<typename T>
T rho_value(matrix<T>& A, matrix<T>& x)
{
	T rho = 0;
	auto temp = A * x;
	for (int i = 0; i < x.rows; ++i)
	{
		rho += temp(i, 0) * x(i, 0);
	}
	return rho;
}

template<typename T>
matrix<T> b_vector(matrix<T>& A)
{
	auto b = matrix<T>(1, A.cols);
	for (int i = 0; i < A.cols; ++i)
	{
		b[i] = A(A.rows - 1, i);
	}
	return b;
}

//int main()
//{
//	using T = double;
//	T epsilon = 1e-8;
//	matrix<T> A = matrix<T>(4, 4, { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 });
//	matrix<T> B = matrix<T>(4, 4, {1.50,0,-0.43,-0.75,0,3,0.87,-0.5,-0.43,0.87,2.9,-0.22,-0.75,-0.5,-0.22,2.6 });
//	std::cout << B << std::endl;
//	//auto Am = A.minor(2);
//	//std::cout << Am << std::endl;
//	//Hessenberg_form(B);
//	//std::cout << B << std::endl;
//	//auto QR = QR_matrices(B);
//	//auto Q = QR.first;
//	//auto R = QR.second;
//	//std::cout << "Q: \n" << Q << '\n' << "R: \n" << R << std::endl;
//	//auto B1 = Q.tranpose() * R;
//	//std::cout << "B1: \n" << B1 << std::endl;
//	auto C = matrix<T>(2, 2, { 1, 6, 1, 2 });
//	QR_eigenvalues(B, epsilon);
//	QR_eigenvalues(B, epsilon, false);
//	QR_eigenvalues(B, epsilon, true, false);
//	auto eigenvalues = QR_eigenvalues(B,epsilon, false, false);
//	std::cout << eigenvalues << std::endl;
//	//auto eigenvectors = reverse_iter_eigenvectors(B, eigenvalues, epsilon);
//	//std::cout << eigenvectors << std::endl;
//	//auto eigensystem = Rayleigh(C, epsilon);
//	//std::cout << eigensystem.first << std::endl;
//	//std::cout << eigensystem.second << std::endl;
//	//auto eigenvalues1 = QR_eigenvalues(C, epsilon);
//	//std::cout << eigenvalues1 << std::endl;
//	//auto eigenvectors1 = reverse_iter_eigenvectors(C, eigenvalues1, epsilon);
//	//std::cout << eigenvectors1 << std::endl;
//}
