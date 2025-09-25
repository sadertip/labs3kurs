#pragma once
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include "Matrices_template.h"

template<typename T>
matrix<T>::~matrix()
{
    delete[] data;
    data = nullptr;
    rows = 0;
    cols = 0;
}

template<typename T>
matrix<T>::matrix(size_t m, size_t n, T value)
    :matrix<T>(m, n)
{
    for (int i = 0; i < m * n; ++i)
    {
        data[i] = value;
    }
}

template<typename T>
matrix<T>::matrix(const matrix<T>& M)
    :matrix<T>(M.row_size(), M.col_size())
{
    auto m = M.row_size();
    auto n = M.col_size();
    for (int i = 0; i < m * n; ++i)
    {
        data[i] = M[i];
    }
}

template<typename T>
matrix<T>::matrix(size_t m, size_t n, std::initializer_list<T> initer)
    :matrix<T>(m, n)
{
    _ASSERT(m * n <= initer.size());
    int i = 0;
    for (auto x : initer)
    {
        data[i] = x;
        ++i;
    }
}

template<typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& M)
{
    auto new_matrix = matrix<T>(M);
    this->swap(new_matrix);
    return *this;
}

template<typename T>
const T& matrix<T>::operator[](int pos) const
{
    _ASSERT(pos < rows * cols);
    return *(data + pos);
}

template<typename T>
T& matrix<T>::operator[](int pos)
{
    _ASSERT(pos < rows * cols);
    return *(data + pos);
}

template<typename T>
matrix<T> matrix<T>::operator*(const matrix<T>& A) const&
{
    _ASSERT(this->cols == A.rows);
    matrix<T> res = matrix<T>(this->rows, A.cols, 0);
    for (int i = 0; i < this->rows; ++i)
        for (int k = 0; k < this->cols; ++k)
            for (int j = 0; j < A.cols; ++j)
                res[res.pos(i, j)] += (*this)[(*this).pos(i, k)] * A[A.pos(k, j)];
    return res;
}

template<typename T>
matrix<T> matrix<T>::row_mult(int m_row, T k) &
{
    for (int j = m_row * (this->cols); j < (m_row + 1) * (this->cols); ++j)
    {
        (*this)[j] *= k;
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_linsum(int m_row, const int n_row, T mult) &
{
    for (int j = 0; j < this->cols; ++j)
    {
        (*this)[m_row * (this->cols) + j] += mult * (*this)[n_row * (this->cols) + j];
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_swap(int m_row, int n_row) &
{
    T* temp_row = new T[this->cols];
    for (int j = 0; j < this->cols; ++j)
    {
        temp_row[j] = (*this)[m_row * this->cols + j];
    }
    for (int j = 0; j < this->cols; ++j)
    {
        (*this)[m_row * this->cols + j] = (*this)[n_row * this->cols + j];
        (*this)[n_row * this->cols + j] = temp_row[j];
    }
    delete[] temp_row;
    return *this;
}

template<typename T>
void matrix<T>::swap(matrix<T>& M)
{
    std::swap(this->rows, M.rows);
    std::swap(this->cols, M.cols);
    std::swap(this->data, M.data);
}

template<typename T>
int matrix<T>::pos(int i, int j) const &
{
    return i * this->cols + j;
}

template<typename T>
T& matrix<T>::operator()(size_t i, size_t j) &
{
    return data[i * cols + j];
}

template<typename T>
const T& matrix<T>::operator()(size_t i, size_t j) const &
{
    return data[i * cols + j];
}

template<typename T>
size_t matrix<T>::row_size() const &
{
    return this->rows;
}

template<typename T>
size_t matrix<T>::col_size() const &
{
    return this->cols;
}

template<typename T>
T* matrix<T>::begin()
{
    return data;
}

template<typename T>
T* matrix<T>::end()
{
    return data + rows * cols;
}

template<typename T>
void matrix<T>::readFromFile(std::ifstream& file)
{
    if (!(file >> rows >> cols)) {
        throw std::runtime_error("Matrix size read error");
    }
    delete[] data;
    data = new T[rows * cols];

    for (size_t i = 0; i < rows * cols; ++i) {
        if (!(file >> data[i])) {
            delete[] data;
            data = nullptr;
            rows = cols = 0;
            throw std::runtime_error("Matrix data read error");
        }
    }
}

template<typename T>
void matrix<T>::writeToFile(const std::string& filename)
{
    std::ofstream file(filename);
    file << rows << ' ' << cols << "\n";

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            file << data[i * cols + j] << " ";
        }
        file << "\n";
    }

    file.close();
}

template<typename T>
matrix<T>::matrix(const std::string& filename)
    :matrix()
{
    std::ifstream file(filename);
    try
    {
        (*this).readFromFile(file);

        file.close(); // Закрываем файл после чтения

    }
    catch (const std::exception& e) {
        //std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "File opening error" << std::endl;
        if (file.is_open())
        {
            file.close();
        }
        //return 1;
    }
    
}

template<typename T>
std::ostream& operator<<(std::ostream& out, matrix<T>& M)
{
    auto m = M.row_size();
    auto n = M.col_size();
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            out << M[i * n + j] << ' ';
        }
        out << '\n';
    }
    return out;
}

template<typename T>
matrix<T> eye(size_t m, T value)
{
    matrix<T> E = matrix<T>(m, m, 0);
    if (value != 0)
        for (int i = 0; i < m; ++i)
            E(i ,i) = value;
    return E;
}

//int main()
//{
//    using T = double;
//    matrix<T> A = matrix<T>(2, 3, { 1, 2, 3, 4, 5, 6});
//    matrix<T> B = matrix<T>(3, 2, { 9, 8, 7, 6, 5,4 });
//    std::cout << A << '\n' << B << std::endl;
//    auto C = A * B;
//    std::cout << C << std::endl;
//    auto E = eye<T>(2, 1);
//    auto temp = C * E;
//    C = E * C;
//    std::cout << temp << '\n' << C << std::endl;
//    return 0;
//}