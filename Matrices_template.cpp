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
matrix<T>::matrix(matrix&& M) noexcept
{
    this->swap(M);
}

template<typename T>
matrix<T>& matrix<T>::operator=(const matrix<T>& M)
{
    auto new_matrix = matrix<T>(M);
    this->swap(new_matrix);
    return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator=(matrix<T>&& M) noexcept
{
    this->swap(M);
    return (*this);
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
matrix<T> matrix<T>::operator*(T k) const {
    matrix<T> result = matrix<T>(this->rows, this->cols);
    for (size_t i = 0; i < this->rows * this->cols; ++i) {
        result.data[i] = data[i] * k;
    }
    return result;
}
template<typename T>
matrix<T>matrix<T>::operator+(const matrix<T>& M) const& {
    if (this->rows != M.rows || this->cols != M.cols) {
        throw std::runtime_error("Matrix dimensions must agree for addition.");
    }
    matrix<T> result = matrix<T>(this->rows, this->cols);
    for (size_t i = 0; i < this->rows * this->cols; ++i) {
        result.data[i] = data[i] + M.data[i];
    }
    return result;
}

// Вычитание матриц/векторов
template<typename T>
matrix<T> matrix<T>::operator-(const matrix<T>& M) const& {
    if (rows != M.rows || cols != M.cols) {
        throw std::runtime_error("Matrix dimensions must agree for subtraction.");
    }
    matrix<T> result = matrix<T>(this->rows, this->cols);
    for (size_t i = 0; i < this->rows * this->cols; ++i) {
        result.data[i] = data[i] - M.data[i];
    }
    return result;
}
template<typename T>
matrix<T> matrix<T>::row_mult(int m_row, T k)&
{
    for (int j = m_row * (this->cols); j < (m_row + 1) * (this->cols); ++j)
    {
        (*this)[j] *= k;
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_linsum(int m_row, const int n_row, T mult)&
{
    for (int j = 0; j < this->cols; ++j)
    {
        (*this)[m_row * (this->cols) + j] += mult * (*this)[n_row * (this->cols) + j];
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_swap(int m_row, int n_row)&
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
matrix<T> matrix<T>::tranpose()&
{
    matrix<T> AT = matrix(this->cols, this->rows);
    for (int j = 0; j < this->cols; ++j)
    {
        for (int i = 0; i < this->rows; ++i)
            AT(j, i) = (*this)(i, j);
    }
    return AT;
}

template<typename T>
matrix<T> matrix<T>::minor(int n)
{
    matrix<T> minor = matrix<T>(n, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            minor(i, j) = (*this)(i, j);
    return minor;
}

template<typename T>
void matrix<T>::swap(matrix<T>& M)
{
    std::swap(this->rows, M.rows);
    std::swap(this->cols, M.cols);
    std::swap(this->data, M.data);
}

template<typename T>
int matrix<T>::pos(int i, int j) const&
{
    return i * this->cols + j;
}

template<typename T>
T& matrix<T>::operator()(size_t i, size_t j)&
{
    return data[i * cols + j];
}

template<typename T>
const T& matrix<T>::operator()(size_t i, size_t j) const&
{
    return data[i * cols + j];
}

template<typename T>
size_t matrix<T>::row_size() const&
{
    return this->rows;
}

template<typename T>
size_t matrix<T>::col_size() const&
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
        std::cerr << "Error: " << e.what() << std::endl;
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
            out << M(i,j) << ' ';
        }
        out << '\n';
    }
    return out;
}
template<typename T>
std::ostream& operator<<(std::ostream& out, const matrix<T>& M)
{
    auto m = M.row_size();
    auto n = M.col_size();
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            out << M(i, j) << ' ';
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
            E(i, i) = value;
    return E;
}
template<typename T>
matrix<T> matrix<T>::get_D() const {

    if (rows != cols) {
        throw std::runtime_error("Matrix must be square for decomposition (D).");
    }
    matrix<T> D(rows, cols, 0.0);

    for (size_t i = 0; i < rows; ++i) {

        D.data[i * cols + i] = data[i * cols + i];
    }
    return D;
}

template<typename T>
matrix<T> matrix<T>::get_L() const {

    if (rows != cols) {
        throw std::runtime_error("Matrix must be square for decomposition (L).");
    }
    matrix<T> L(rows, cols, 0.0);

    for (size_t i = 1; i < rows; ++i) {
        for (size_t j = 0; j < i; ++j) {

            L.data[i * cols + j] = data[i * cols + j];
        }
    }
    return L;
}

template<typename T>
matrix<T> matrix<T>::get_U() const {

    if (rows != cols) {
        throw std::runtime_error("Matrix must be square for decomposition (U).");
    }
    matrix<T> U(rows, cols, 0.0);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = i + 1; j < cols; ++j) {
            U.data[i * cols + j] = data[i * cols + j];
        }
    }
    return U;
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