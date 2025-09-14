#pragma once
class matrix
{
	size_t rows = 0;
	size_t cols = 0;
	double* data = nullptr;
public:
	matrix() : rows(0), cols(0), data(nullptr) {}
	matrix(size_t m, size_t n) : rows(m), cols(n), data(new double[m * n]) {}
	matrix(size_t m, size_t n, double value);
	matrix(const matrix& A);
	matrix(size_t m, size_t n, std::initializer_list<double> arr);
	const double& operator[](int position) const;
	double& operator[](int position);
	matrix row_mult(int m, double k);

	size_t row_size() const;
	size_t col_size() const;

	~matrix();

	double* begin();

	double* end();

	void readFromFile(const std::string& filename);

	void writeToFile(const std::string& filename);
};

std::ostream& operator<<(std::ostream& out, matrix& A)
{
	auto m = A.row_size();
	auto n = A.col_size();
	for(int i = 0; i < m; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			out << A[i * n + j] << ' ';
		}
		out << '\n';
	}
	return out;
}