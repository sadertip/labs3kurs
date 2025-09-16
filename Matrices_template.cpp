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
    :matrix(m, n)
{
    for (int i = 0; i < m * n; ++i)
    {
        data[i] = value;
    }
}

template<typename T>
matrix<T>::matrix(const matrix& M)
    :matrix(M.row_size(), M.col_size())
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
    :matrix(m, n)
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
matrix<T>& matrix<T>::operator=(const matrix& M)
{
    auto new_matrix = matrix(M);
    this->swap(new_matrix);
    return *this;
}

template<typename T>
const T& matrix<T>::operator[](int position) const
{
    _ASSERT(position < rows * cols);
    return *(data + position);
}

template<typename T>
T& matrix<T>::operator[](int position)
{
    _ASSERT(position < rows * cols);
    return *(data + position);
}

template<typename T>
matrix<T> matrix<T>::row_mult(int m_row, T k)
{
    for (int j = m_row * (this->cols); j < (m_row + 1) * (this->cols); ++j)
    {
        (*this)[j] *= k;
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_linsum(int m_row, const int n_row, T mult)
{
    for (int j = 0; j < this->cols; ++j)
    {
        (*this)[m_row * (this->cols) + j] += mult * (*this)[n_row * (this->cols) + j];
    }
    return *this;
}

template<typename T>
matrix<T> matrix<T>::row_swap(int m_row, int n_row)
{
    T* temp_row = new T[this->cols];
    for (int j = 0; j < this->cols; ++j)
    {
        *(temp_row + j) = (*this)[m_row * this->cols + j];
    }
    for (int j = 0; j < this->cols; ++j)
    {
        (*this)[m_row * this->cols + j] = (*this)[n_row * this->cols + j];
        (*this)[n_row * this->cols + j] = *(temp_row + j);
    }
    return *this;
}

template<typename T>
void matrix<T>::swap(matrix& M)
{
    std::swap(this->rows, M.rows);
    std::swap(this->cols, M.cols);
    std::swap(this->data, M.data);
}

template<typename T>
size_t matrix<T>::row_size() const
{
    return this->rows;
}

template<typename T>
size_t matrix<T>::col_size() const
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
void matrix<T>::readFromFile(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Не удалось открыть файл: " + filename);
    }

    // Читаем размеры матрицы
    file >> rows >> cols;

    if (rows == 0 || cols == 0) {
        throw std::runtime_error("Некорректные размеры матрицы");
    }

    // Освобождаем старые данные
    delete[] data;
    data = new T[rows * cols];

    // Читаем элементы матрицы
    for (size_t i = 0; i < rows * cols; ++i) {
        if (!(file >> data[i])) {
            delete[] data;
            data = nullptr;
            rows = cols = 0;
            throw std::runtime_error("Ошибка чтения данных из файла");
        }
    }

    file.close();
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
void gauss_step(matrix<T>& M, int row)
{
    auto first_element = M.data[row * M.cols];
    auto main_elem = first_element;
    for (size_t i = first_element; i <= row * M.cols; ++i)
    {
        M.data[i] /= main_elem;
    }
    for (size_t i = row; i <= M.rows; ++i)
    {
        M[i] = M[i] - M[row];
    }

}

int main()
{
    using T = double;
    matrix<T> A = matrix<T>(2, 2, 1);
    matrix<T> B = matrix<T>(3, 2, { 1, 2, 3, 4, 5, 6 });
    std::cout << B << std::endl;
    B.row_mult(0, 3);
    A = B;
    A.row_linsum(2, 0, 2);
    std::cout << B << '\n' << A << std::endl;
    matrix<T> C = matrix<T>(4, 4, 1);
    std::cout << C << std::endl;
    C.row_mult(1, 2);
    C.row_mult(2, 3);
    C.row_linsum(3, 1, 2);
    C.row_linsum(3, 2, 1);
    std::cout << C << std::endl;
}

