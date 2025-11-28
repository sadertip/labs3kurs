#pragma once
//#include "Nets.cpp"
#include "Matrices_template.h"

template<typename T>
class Lagrange_polynom
{
public:
	matrix<T> points;
	matrix<T> divided_differences;
	size_t degree;
	Lagrange_polynom() : points(), divided_differences(), degree(0) {}
	Lagrange_polynom(matrix<T>& grid, T(*f)(T));
	Lagrange_polynom(matrix<T>& grid, matrix<T>& values);

	T operator()(T x)&;


	~Lagrange_polynom();
};

template<typename T>
class Spline
{
public:
	matrix<T> points;
	matrix<T> a_Coeffs;
	matrix<T> b_Coeffs;
	matrix<T> c_Coeffs;
	matrix<T> d_Coeffs;
	Spline() : points(), a_Coeffs(), b_Coeffs(), c_Coeffs(), d_Coeffs() {}
	Spline(matrix<T>& grid, T(*f)(T));
	Spline(matrix<T>& grid, matrix<T>& values);

	T operator()(T x)&;


	~Spline();
};

template<typename T>
matrix<T> generate_values(matrix<T>& grid, T(*f)(T));

template<typename T>
matrix<T> generate_values(matrix<T>& grid, Lagrange_polynom<T>& L);

template<typename T>
matrix<T> generate_values(matrix<T>& grid, Spline<T>& S);

template<typename T>
matrix<T> tridiag_solve(int n, matrix<T>& a, matrix<T>& b, matrix<T>& c, matrix<T>& d);
//template<typename T>
//void write_to_file(matrix<T>& values, )
template<typename T>
T func1(T x);

template<typename T>
T func2(T x);

template<typename T>
T func3(T x);

template<typename T>
T func4(T x);

template<typename T>
T Const(T x);

template<typename T>
T Lin(T x);

template<typename T>
T Runge(T x);