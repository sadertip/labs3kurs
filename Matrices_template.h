#pragma once
#include <iostream>
#include <fstream>

template<typename T>
class matrix
{
public:
	size_t rows = 0;
	size_t cols = 0;
	T* data = nullptr;
	matrix() : rows(0), cols(0), data(nullptr) {}
	matrix(size_t m, size_t n) : rows(m), cols(n), data(new T[m * n]) {}
	matrix(size_t m, size_t n, T value);
	matrix(const matrix<T>& M);
	matrix(size_t m, size_t n, std::initializer_list<T> arr);
	matrix(matrix&& M) noexcept;

	matrix<T>& operator=(const matrix<T>& M);
	matrix<T>& operator=(matrix<T>&& M) noexcept;

	const T& operator[](int pos) const;
	T& operator[](int pos);

	matrix<T> operator*(const matrix<T>& A) const&;
	//matrix<T> operator*(const matrix<T>& A)&&;

	matrix<T> row_mult(int m_row, T k) &;
	matrix<T> row_linsum(int m_row, const int n_row, T mult) &;
	matrix<T> row_swap(int m_row, int n_row) &;

	void swap(matrix<T>& M);

	inline int pos(int i, int j) const &;

	T& operator()(size_t i, size_t j) &;

	const T& operator()(size_t i, size_t j) const &;

	size_t row_size() const &;
	size_t col_size() const &;

	~matrix();

	T* begin();

	T* end();

	void readFromFile(std::ifstream& file);

	void writeToFile(const std::string& filename);

	matrix(const std::string& filename);
};

template<typename T>
std::ostream& operator<<(std::ostream& out, matrix<T>& M);

template<typename T>
matrix<T> eye(size_t m, T value);
