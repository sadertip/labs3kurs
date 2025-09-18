#pragma once

#include "Matrices_template.h"

matrix<double> QR_decomposition(matrix<double>& A, matrix<double>& b);

//double C_ij(matrix<double>& A, int i, int j);
//double S_ij(matrix<double>& A, int i, int j);

std::pair<double, double> Coeffs_ij(matrix<double>& A, int i, int j);

void QR_row_step(matrix<double>& A, matrix<double>& b, double c_ij, double s_ij, int row, int step);
void QR_col_step(matrix<double>& A, matrix<double>& b, matrix<double>& Q, int j);

