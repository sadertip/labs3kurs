#pragma once
#include <iostream>
#include "Matrices_template.cpp"
#include "Gaus.h"

template<typename T>
void pivot(matrix<T>& M, matrix<T>& b, int row)
{
    auto maxi = M[row * M.cols + row];
    size_t maxi_elem_num = row * M.cols + row;
    for (int i = row * M.cols + row; i <= M.cols * M.cols - 1; i += M.cols)
    {
        if (maxi < M[i])
        {
            maxi = M[i];
            maxi_elem_num = i;
        }
    }
    if (maxi_elem_num != row)
    {
        M.row_swap(row, maxi_elem_num / M.cols);
        b.row_swap(row, maxi_elem_num / M.cols);
    }
}

template<typename T>
void gauss_step(matrix<T>& M, matrix<T>& b, int row)
{
    pivot(M, b, row);
    std::cout << M << std::endl;
    std::cout << b << std::endl;
    size_t main_elem_num = row * M.cols;
    auto main_elem = M.data[row * M.cols];

    for (auto i = row * M.cols; i <= row * M.cols + M.cols - 1; ++i)
    {
        if (M[i] != 0)
        {
            main_elem = M[i];
            main_elem_num = i;
            break;
        }
    }
    M.row_mult(row, 1 / main_elem);
    b.row_mult(row, 1 / main_elem);
    for (int i = row + 1; i < M.rows; ++i)
    {
        auto k = -M[main_elem_num - M.cols * row + M.cols * i];
        M.row_linsum(i, row, k);
        b.row_linsum(i, row, k);
    }

}

template<typename T>
matrix<T> gauss_back(matrix<T>& A, matrix<T>& b)
{
    size_t n = A.rows;
    matrix<T> x(n, 1);

    for (int i = n - 1; i >= 0; --i)
    {
        T sum = 0.0;
        for (size_t j = i + 1; j < n; ++j)
        {
            sum += A[A.cols * i + j] * x[j];
        }

        if (std::abs(A[A.cols * i + i]) < 1e-9)
        {
            throw std::runtime_error("Division by zero occurred. matrix<T> may not be invertible.");
        }

        x[i] = (b[i] - sum) / A[A.cols * i + i];
    }

    return x;
}

template<typename T>
matrix<T> Gauss(matrix<T>& A, matrix<T>& b)
{

    for (size_t i = 0; i < A.rows - 1; ++i)
    {
        gauss_step(A, b, i);
    }
    matrix<T> x = gauss_back(A, b);
    return x;
}

//int main()
//{
//    using T = double;
//    std::ifstream file("input.txt");
//    matrix<T> A, b;
//    try {
//
//
//        // Читаем первую матрицу (A)
//        A.readFromFile(file);
//
//        // Читаем вторую матрицу (b)
//        b.readFromFile(file);
//
//        file.close(); // Закрываем файл после чтения
//
//    }
//    catch (const std::exception& e) {
//        std::cerr << "Ошибка: " << e.what() << std::endl;
//        if (file.is_open()) {
//            file.close();
//        }
//        return 1;
//    }
//    std::cout << A << std::endl;
//    std::cout << b << std::endl;
//    matrix<T> x = Gauss(A, b);
//    std::cout << x << std::endl;
//}