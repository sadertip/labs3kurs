#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include "Matrices.h"


matrix::~matrix()
{
    delete[] data;
    data = nullptr;
    rows = 0;
    cols = 0;
}

matrix::matrix(size_t m, size_t n, double value)
    :matrix(m, n)
{
    for (int i = 0; i < m * n; ++i)
    {
        data[i] = value;
    }
}

matrix::matrix(const matrix& A)
    :matrix(A.row_size(), A.col_size())
{
    auto m = A.row_size();
    auto n = A.col_size();
    for (int i = 0; i < m * n; ++i)
    {
        data[i] = A[i];
    }
}

matrix::matrix(size_t m, size_t n, std::initializer_list<double> initer)
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

const double& matrix::operator[](int position) const
{
    _ASSERT(position < rows * cols);
    return *(data + position);
}

double& matrix::operator[](int position)
{
    _ASSERT(position < rows * cols);
    return *(data + position);
}

matrix matrix::row_mult(int m, double k)
{
    for (int j = m * this->cols; j < m + this->cols; ++j)
    {
        (*this)[j] *= k;
    }
    return *this;
}

size_t matrix::row_size() const
{
    return this->rows;
}

size_t matrix::col_size() const
{
    return this->cols;
}

double* matrix::begin()
{
    return data;
}

double* matrix::end()
{
    return data + rows * cols;
}
void matrix::readFromFile(const std::string& filename)
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
    data = new double[rows * cols];

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
void matrix::writeToFile(const std::string& filename)
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
void gauss_step(matrix& A, int row)
{
    auto first_element = A.data[row * A.cols];
    auto main_elem = first_element;
    for (size_t i = first_element; i <= row * A.cols; ++i)
    {
        A.data[i] /= main_elem;
    }
    for (size_t i = row; i <= A.rows; ++i)
    {
        A[i] = A[i] - A[row];
    }

}

int main()
{
    matrix A = matrix(2, 2, 1);
    matrix B = matrix(2, 3, { 1, 2, 3, 4, 5, 6 });
    B.row_mult(2, 2);
    std::cout << B << std::endl;
}

