#pragma once
class matrix
{
public:
	size_t rows = 0;
	size_t cols = 0;
	double* data = nullptr;
	matrix() : rows(0), cols(0), data(nullptr) {}
	matrix(size_t m, size_t n) : rows(m), cols(n), data(new double[m * n]) {}
	matrix(size_t m, size_t n, double value);
	matrix(const matrix& M);
	matrix(size_t m, size_t n, std::initializer_list<double> arr);
	
	matrix& operator=(const matrix& M);
	
	const double& operator[](int position) const;
	double& operator[](int position);

	matrix row_mult(int m_row, double k);
	matrix row_linsum(int m_row, const int n_row, double mult);
	matrix row_swap(int m_row, int n_row);
	
	void swap(matrix& M);

	size_t row_size() const;
	size_t col_size() const;

	~matrix();

	double* begin();

	double* end();

	void readFromFile(const std::string& filename);

	void writeToFile(const std::string& filename);
};

std::ostream& operator<<(std::ostream& out, matrix& M);