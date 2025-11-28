#pragma once
#include <iostream>
#include "Interpolation.h"
#include "Nets.cpp"
//#include "Matrices_template.cpp"
#include <cmath>

#define SQR(x) (x*x)
#define CUBE(x) (x*x*x)
#define EPSILON 10e-6

template<typename T>
Lagrange_polynom<T>::Lagrange_polynom(matrix<T>& grid, T(*f)(T))
{
    //degree = grid.cols;
    //points = matrix<T>(grid);
    //divided_differences(grid.rows, grid.cols);

    auto prev_col = generate_values<T>(grid, f);
    (*this) = Lagrange_polynom(grid, prev_col);
    //auto cur_col = matrix<T>(grid.rows, grid.cols - 1);
   
    //divided_differences[0] = f(grid[0]);
    //for (int k = 1; k < degree; ++k)
    //{
    //    for (int i = 0; i < degree - k; ++i)
    //    {
    //        cur_col[i] = (prev_col[i + 1] - prev_col[i]) / (grid[i + k] - grid[i]);
    //    }
    //    divided_differences[k] = cur_col[0];
    //    prev_col.cols -= 1;
    //    prev_col = cur_col;
    //}
}

template<typename T>
Lagrange_polynom<T>::Lagrange_polynom(matrix<T>& grid, matrix<T>& values)
{
    degree = grid.cols;
    points = matrix<T>(grid);
    divided_differences = matrix<T>(grid.rows, grid.cols);

    auto prev_col = matrix<T>(values);
    auto cur_col = matrix<T>(grid.rows, grid.cols - 1);
    //std::cout << "cur_col:\n" << cur_col << std::endl;

    divided_differences[0] = values[0];
    for (int k = 1; k < degree; ++k)
    {
        for (int i = 0; i < degree - k; ++i)
        {
            cur_col[i] = (prev_col[i + 1] - prev_col[i]) / (grid[i + k] - grid[i]);
        }
        //std::cout << k <<" cur_col:\n" << cur_col << std::endl;

        divided_differences[k] = cur_col[0];
        //std::cout << k << " div_diff:\n" << divided_differences << std::endl;
        prev_col.cols -= 1;
        prev_col = cur_col;
        cur_col.cols -= 1;
    }
}

template<typename T>
T Lagrange_polynom<T>::operator()(T x)&
{
    T ans = 0;
    for (int i = 0; i < degree - 1; ++i)
    {
        T prod = this->divided_differences[i];
        for (int k = 0; k < i; ++k)
        {
            prod *= x - this->points[k];
        }
        ans += prod;
    }
    return ans;
}

template<typename T>
Lagrange_polynom<T>::~Lagrange_polynom()
{
    this->points.~matrix();
    this->divided_differences.~matrix();
}

template<typename T>
Spline<T>::Spline(matrix<T>& grid, T(*f)(T))
{
    auto prev_col = generate_values<T>(grid, f);
    (*this) = Spline(grid, prev_col);
}

template<typename T>
Spline<T>::Spline(matrix<T>& grid, matrix<T>& values)
{
    points = matrix<T>(grid);
    int n = grid.cols - 1;

    auto h = matrix<T>(1, n);
    auto g = matrix<T>(1, n);

    for (int i = 0; i < n; ++i)
        h[i] = grid[i + 1] - grid[i];
    //std::cout << "h:\n" << h << std::endl;
    for (int i = 0; i < n; ++i)
        g[i] = (values[i + 1] - values[i]) / h[i];
    //std::cout << "g:\n" << g << std::endl;

    a_Coeffs = matrix<T>(1, n);
    b_Coeffs = matrix<T>(1, n);
    c_Coeffs = matrix<T>(1, n);
    d_Coeffs = matrix<T>(1, n);

    auto a = matrix<T>(1, n);
    auto b = matrix<T>(1, n);
    auto c = matrix<T>(1, n);
    auto d = matrix<T>(1, n);

    a[0] = 0;
    for (int i = 1; i < n; ++i)
        a[i] = h[i];
    //std::cout << "a:\n" << a << std::endl;
    b[0] = 1;
    for (int i = 1; i < n; ++i)
        b[i] = -2*(h[i - 1] + h[i]);
    //std::cout << "b:\n" << b << std::endl;
    c[0] = 0;
    c[n - 1] = 0;
    for (int i = 1; i < n-1; ++i)
        c[i] = h[i];
    //std::cout << "c:\n" << c << std::endl;
    d[0] = 0;
    for (int i = 1; i < n ; ++i)
        d[i] = -3*(g[i] - g[i-1]);
    //std::cout << "d:\n" << d << std::endl;
    auto x = tridiag_solve(n, a, b, c, d);
    //std::cout << "x:\n" << x << std::endl;
    c_Coeffs = x;
    for (int i = 0; i < n - 1; ++i)
        b_Coeffs[i] = g[i] - (c_Coeffs[i + 1] + 2 * c_Coeffs[i]) * h[i] / 3;
    b_Coeffs[n - 1] = g[n - 1] - 2 * c_Coeffs[n - 1] * h[n - 1] / 3;
    //std::cout << "b_Co:\n" << b_Coeffs << std::endl;
    for (int i = 0; i < n - 1; ++i)
        d_Coeffs[i] = (c_Coeffs[i + 1] - c_Coeffs[i]) / (3*h[i]);
    d_Coeffs[n - 1] = -c_Coeffs[n-1] / (3 * h[n-1]);
    //std::cout << "d_Co:\n" << d_Coeffs << std::endl;
    for (int i = 0; i < n; ++i)
        a_Coeffs[i] = values[i];
    //std::cout << "a_Co:\n" << a_Coeffs << std::endl;
}

template<typename T>
T Spline<T>::operator()(T x)&
{
    for (int i = 0; i < points.cols - 1; ++i)
    {
        //std::cout << points << std::endl;
        //std::cout << points[i] - EPSILON << std::endl;
        //std::cout << points[i+1] + EPSILON << std::endl;
        if ((x >= (points[i] - EPSILON)) && (x <= (points[i + 1] + EPSILON)))
        {
            T ans = 0;
            ans += a_Coeffs[i];
            ans += b_Coeffs[i] * (x - points[i]);
            ans += c_Coeffs[i] * SQR((x - points[i]));
            ans += d_Coeffs[i] * CUBE((x - points[i]));
            return ans;
        }
    }
}

template<typename T>
Spline<T>::~Spline()
{
    this->a_Coeffs.~matrix();
    this->b_Coeffs.~matrix();
    this->c_Coeffs.~matrix();
    this->d_Coeffs.~matrix();
    this->points.~matrix();
}

template<typename T>
matrix<T> generate_values(matrix<T>& grid, T(*f)(T))
{
    auto ans = matrix<T>(grid.rows, grid.cols);
    for (int i = 0; i < grid.cols; ++i)
    {
        ans[i] = f(grid[i]);
    }
    return ans;
}

template<typename T>
matrix<T> generate_values(matrix<T>& grid, Lagrange_polynom<T>& L)
{
    auto ans = matrix<T>(grid.rows, grid.cols);
    for (int i = 0; i < grid.cols; ++i)
    {
        ans[i] = L(grid[i]);
    }
    return ans;
}

template<typename T>
matrix<T> generate_values(matrix<T>& grid, Spline<T>& S)
{
    auto ans = matrix<T>(grid.rows, grid.cols);
    for (int i = 0; i < grid.cols; ++i)
    {
        ans[i] = S(grid[i]);
    }
    return ans;
}

template<typename T>
matrix<T> tridiag_solve(int n, matrix<T>& a, matrix<T>& b, matrix<T>& c, matrix<T>& d)
{
    auto x = matrix<T>(1, n);

    auto alpha = matrix<T>(1, n+1);
    alpha[0] = 0;

    auto beta = matrix<T>(1, n+1);
    beta[0] = 0;

    for (int i = 0; i < n; ++i)
    {
        T den = b[i] - a[i] * alpha[i];
        alpha[i + 1] = c[i] / den;
        //std::cout << alpha[i+1] << std::endl;
        beta[i + 1] = (d[i] + a[i] * beta[i]) / den;
        //std::cout << beta[i + 1] << std::endl;
    }
    //std::cout << alpha << std::endl;
    //std::cout << beta << std::endl;
    x[n - 1] = beta[n];
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
    }
    return x;
}

template<typename T>
T func1(T x)
{
    return SQR(x);
}

template<typename T>
T func2(T x)
{
    return 1 / (1 + SQR(x));
}

template<typename T>
T func3(T x)
{
    return 1 / (std::atan(1+ 10*SQR(x)));
}

template<typename T>
T func4(T x)
{
    T a = pow((4 * CUBE(x) + 2 * SQR(x) - 4 * x + 2), sqrt(2));
    T b = std::asin(1 / (5 + x - SQR(x)));
    return a + b - 5;
}

template<typename T>
T Const(T x)
{
    return 1;
}

template<typename T>
T Lin(T x)
{
    return x;
}

template<typename T>
T Runge(T x)
{
    return 1 / (1 + 25 * SQR(x));
}
