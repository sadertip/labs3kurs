#pragma once
#include "Matrices_template.h"

template<typename T>
matrix<T> generate_uniform_grid_matrix(T a, T b, int n);

template<typename T>
matrix<T> generate_chebyshev_grid_matrix(T a, T b, int n);