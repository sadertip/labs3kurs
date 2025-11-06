#pragma once
#include "Matrices_template.h"

template<typename T>
struct SorComponents {
    matrix<T> G1;
    matrix<T> G2;
    matrix<T> C;
    matrix<T> CU;
    matrix<T> CL;
    SorComponents(const matrix<T>& g1, const matrix<T>& g2, const matrix<T>& c,const matrix<T>& cu, const matrix<T>& cl)
        : G1(g1), G2(g2), C(c), CU(cu), CL(cl){}
};
template <typename T>
T sim_rho(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T tau);
template<typename T>
int est_iter(T q, T epsilon, T rho_0);
template<typename T>
std::pair<matrix<T>, matrix<T>> get_simple_iteration_components(const matrix<T>& A, const matrix<T>& b, T tau);
template<typename T>
std::pair<matrix<T>, matrix<T>> get_jacobi_components(const matrix<T>& A, const matrix<T>& b);
template <typename T>
matrix<T> simple_iteration_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T tau, T q, int max_iter, T epsilon);
template<typename T>
matrix<T> jacobi_method_elementwise(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T q, int max_iter, T epsilon);
template<typename T>
matrix<T> sor_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T omega, int max_iter, T epsilon);