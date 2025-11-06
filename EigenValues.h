#pragma once
#include "Matrices_template.cpp"

template<typename T>
void Hessenberg_form(matrix<T>& A);

template<typename T>
std::pair<T, T> Alpha_beta(matrix<T>& A, int k, int l);

template<typename T>
void rotation(matrix<T>& A, int k, int l);

template<typename T>
T sigma_value(matrix<T>& A);

template<typename T>
void shift(matrix<T>& A, T sigma);

template<typename T>
void backshift(matrix<T>& A, T sigma);

template<typename T>
std::pair<matrix<T>, matrix<T>> QR_matrices(matrix<T>& A, bool Hessenberg = false);

template<typename T>
void QR_matrices_row_step(matrix<T>& A, T c_ij, T s_ij, int row, int step, bool Hessenberg = false);

template<typename T>
void QR_matrices_col_step(matrix<T>& A, matrix<T>& Q, int j, bool Hessenberg = false);

template<typename T>
matrix<T> QR_eigenvalues(matrix<T>& A, T epsilon, bool Hessenberg = true, bool shifting = true);

template<typename T>
matrix<T> reverse_iter_eigenvectors(matrix<T>& A, matrix<T>& eigenvalues, T epsilon);

template<typename T>
std::pair<matrix<T>, matrix<T>> Rayleigh(matrix<T>& A, T epsilon);

template<typename T>
T rho_value(matrix<T>& A, matrix<T>& x);

template<typename T>
matrix<T> b_vector(matrix<T>& A);