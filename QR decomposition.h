#pragma once

#include "Matrices_template.h"

template<typename T>
matrix<T> QR_decomposition(matrix<T>& A, matrix<T>& b);

//T C_ij(matrix<T>& A, int i, int j);
//T S_ij(matrix<T>& A, int i, int j);

template<typename T>
std::pair<T, T> Coeffs_ij(matrix<T>& A, int i, int j);

template<typename T>
void QR_row_step(matrix<T>& A, matrix<T>& b, T c_ij, T s_ij, int row, int step);

template<typename T>
void QR_col_step(matrix<T>& A, matrix<T>& b, matrix<T>& Q, int j);

