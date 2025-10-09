#pragma once
#include "Matrices_template.h"

template<typename T>
matrix<T> gauss_back(matrix<T>& A, matrix<T>& b);

template<typename T>
void gauss_step(matrix<T>& M, matrix<T>& b, int row);

template<typename T>
matrix<T> Gauss(matrix<T>& A, matrix<T>& b);

