#pragma once
#include "Matrices_template.h"
template<typename T>
T norm1(matrix<T>& A);
template<typename T>
T norm2(matrix<T>& A);
template<typename T>
T cond(matrix<T>& m);
template<typename T>
matrix<T> reverse(const matrix<T>& m);
template<typename T>
T mismatch(matrix<T>& a, matrix<T>& b);
template<typename T>
T norminf(matrix<T>& A);