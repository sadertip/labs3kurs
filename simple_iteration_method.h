#pragma once
#include "Matrices_template.h"

template<typename T>
struct SorComponents {
    matrix<T> CL;
    matrix<T> CU;
    matrix<T> C;
    SorComponents(const matrix<T>& cl, const matrix<T>& cd,
        const matrix<T>& cu, const matrix<T>& c)
        : CL(cl), CU(cu), C(c) {}
};
template<typename T>
std::pair<matrix<T>, matrix<T>> get_simple_iteration_components(const matrix<T>& A, const matrix<T>& b, T tau);
template<typename T>
std::pair<matrix<T>, matrix<T>> get_jacobi_components(const matrix<T>& A, const matrix<T>& b);
template <typename T>
matrix<T> simple_iteration_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T tau, int max_iter, T epsilon);
template<typename T>
matrix<T> jacobi_method_elementwise(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, int max_iter, T epsilon);
template<typename T>
matrix<T> sor_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T omega, int max_iter, T epsilon);