#pragma once
#include <iostream>
#include <cmath>
#include "Matrices_template.cpp"
#include "simple_iteration_method.h"
#include "other.cpp"



template <typename T>
T sim_rho(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T tau) {
    size_t n = A.rows;
    matrix<T> E = eye(n, 1.0);
    matrix<T> B = E - (A * tau);
    matrix<T> f = b * tau;

    matrix<T> x1 = B * x0 + f;
    return norminf(x1 - x0);
}
template<typename T>
T jacobi_rho(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0)
{
    size_t n = A.rows;


    matrix<T> x_curr = x0;

    matrix<T> x_next(n, 1, 0.0);

    for (size_t i = 0; i < n; ++i)
    {
        T a_ii = A.data[i * n + i];
        if (std::abs(a_ii) < 1e-12) {
            throw std::runtime_error("Zero diagonal element encountered during rho0 calculation.");
        }

        T sum = 0;

 
        for (size_t j = 0; j < n; ++j)
        {
            if (i != j) {
                T a_ij = A.data[i * n + j];
                T x_j_k = x_curr.data[j];
                sum += a_ij * x_j_k;
            }
        }

        T b_i = b.data[i];

        x_next.data[i] = (b_i - sum) / a_ii;
    }

    return norminf(x_next - x_curr);

}
template <typename T>
T sor_rho(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T omega) {
    size_t n = A.rows;
    matrix<T> x1 = x0; 
    matrix<T> x0_copy = x0; 

    for (size_t i = 0; i < n; ++i) {
        T a_ii = A.data[i * n + i];
        T sum_terms = 0;

        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                T a_ij = A.data[i * n + j];
                if (j < i) {
                    sum_terms += a_ij * x1.data[j];
                }
                else {
                    sum_terms += a_ij * x0_copy.data[j];
                }
            }
        }
        T b_i = b.data[i];
        T x_i_gs_approx = (b_i - sum_terms) / a_ii;

        x1.data[i] = (1.0 - omega) * x0_copy.data[i] + omega * x_i_gs_approx;
    }

    return norminf(x1 - x0_copy);
}
template <typename T>
int est_iter(T q, T epsilon, T rho_0) {
    T numerator = (static_cast<T>(1.0) - q) * epsilon;
    T Log_Arg = numerator / rho_0;
    if (Log_Arg >= static_cast<T>(1.0)) {
        return 1;
    }

    T k_exact = std::log(Log_Arg) / std::log(q);

    int k_min = static_cast<int>(std::ceil(k_exact));

    return k_min;
}
template<typename T>
std::pair<matrix<T>, matrix<T>> get_simple_iteration_components(const matrix<T>& A, const matrix<T>& b, T tau)
{
    size_t n = A.rows;

    matrix<T> E = eye(n, 1.0);

    matrix<T> tau_A = A * tau;

    matrix<T> C = E - tau_A;

    matrix<T> y = b * tau;

    return { C, y };
}
template<typename T>
std::pair<matrix<T>, matrix<T>> get_jacobi_components(const matrix<T>& A, const matrix<T>& b)
{

    matrix<T> D = A.get_D();
    std::cout << D;
    matrix<T> D_inv = reverse(D);
    

    matrix<T> L_plus_U = A - D;


    matrix<T> D_inv_L_plus_U = D_inv * L_plus_U;


    matrix<T> C = D_inv_L_plus_U * (static_cast<T>(-1.0));


    matrix<T> y = D_inv * b;

    return { C, y };
}

template<typename T>
SorComponents<T> calculate_sor_matrices(const matrix<T>& A, T omega) {
    size_t n = A.rows;
    if (n != A.cols) {
        throw std::runtime_error("Matrix must be square.");
    }

    matrix<T> D = A.get_D();
    matrix<T> L = A.get_L();
    matrix<T> U = A.get_U();
    
    matrix<T> D_inv = reverse(D);

    matrix<T> G1 = (D_inv * L) * (-omega);

    matrix<T> E = eye(n, 1.0);

    matrix<T> term1 = E * (1.0 - omega);
    matrix<T> term2 = (D_inv * U) * (-omega);

    matrix<T> G2 = term1 + term2;

    matrix<T> omega_L = L * omega;

    matrix<T> DL_term = D + omega_L;

    matrix<T> DL_inv = reverse(DL_term);


    matrix<T> one_minus_omega_D = D * (1.0 - omega);
    matrix<T> minus_omega_U = U * (-omega);

    matrix<T> RHS_term = one_minus_omega_D + minus_omega_U;


    matrix<T> C = DL_inv * RHS_term;


    matrix<T> CU = C.get_U();
    matrix<T> CL = C.get_L();

    return SorComponents<T>{ G1, G2, C,CU,CL };
}
template<typename T>
void calculate_and_sum_norms(const matrix<T>& G1, const matrix<T>& G2) {

    T norm_inf_G1 = norminf(G1);
    T norm_inf_G2 = norminf(G2);

    T sum_norm_inf = norm_inf_G1 + norm_inf_G2;

    std::cout << "--- Infinite Norm (||.||_inf) ---\n";
    std::cout << "  ||CU||_inf = " << norm_inf_G1 << "\n";
    std::cout << "  ||CL||_inf = " << norm_inf_G2 << "\n";
    std::cout << "  ||CU||_inf + ||CL||_inf = " << sum_norm_inf << "\n";

    T norm_1_G1 = norm1(G1);
    T norm_1_G2 = norm1(G2);

    T sum_norm_1 = norm_1_G1 + norm_1_G2;

    std::cout << "--- One Norm (||.||_1) ---\n";
    std::cout << "  ||CU||_1 = " << norm_1_G1 << "\n";
    std::cout << "  ||CL||_1 = " << norm_1_G2 << "\n";
    std::cout << "  ||CU||_1 + ||CL||_1 = " << sum_norm_1 << "\n";
}

template <typename T>
matrix<T> simple_iteration_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T tau, T q, int max_iter, T epsilon)
{

    if (A.rows != A.cols || A.rows != b.rows || b.cols != 1 || A.rows != x0.rows || x0.cols != 1) {
        throw std::runtime_error("Invalid matrix/vector dimensions for simple_iteration_method.");
    }

    size_t n = A.rows;
    matrix<T> E = eye(n, 1.0);


    matrix<T> B = E - (A * tau);

    matrix<T> f = b * tau;
    /*std::cout << b << '\n';
    std::cout << tau << '\n';
    std::cout << f << '\n';
    std::cout << x0 << '\n';*/
    matrix<T> x_curr = x0;
    /*std::cout << x_curr << '\n';*/
    matrix<T> x_next;
    const T epsilon_0 = static_cast<T>(1e-16);
    T C4_Threshold = ((static_cast<T>(1.0) - q) / q) * epsilon;
    matrix<T> Ax0_minus_b = A * x0 - b;
    T Norm_R0 = norminf(Ax0_minus_b);
    std::cout << "Starting Simple Iteration Method (tau = " << tau << ")\n";

    for (int k = 0; k < max_iter; ++k) {
        x_next = B * x_curr + f;
        /*if (k == 0)
        {
            std::cout << x0<<'\n';
            std::cout << x_curr<<'\n';
            std::cout << f << '\n';
            std::cout << x_next << '\n';
            std::cout << B;
        }*/
        matrix<T> residual_diff = x_next - x_curr;
        T diff_norm = norminf(residual_diff);

        matrix<T> R_next_vector = A * x_next - b;
        T R_next_norm = norminf(R_next_vector);

        T x_curr_norm = norminf(x_curr);


        T C1_Value = diff_norm / (x_curr_norm + epsilon_0);
        bool Is_C1_Met = (C1_Value <= epsilon);


        T C2_Value = R_next_norm;
        bool Is_C2_Met = (C2_Value <= epsilon);


        T C3_Value = (Norm_R0 < 1e-16) ? 0.0 : R_next_norm / Norm_R0;
        bool Is_C3_Met = (C3_Value <= epsilon);

        /*T C4_Value = diff_norm; 
        bool Is_C4_Met = (C4_Value <= C4_Threshold);*/


        //std::cout << "Iteration " << k + 1 << ":\n";
        //std::cout << " C1 (Rel. Diff) = " << C1_Value << "\n";
        //std::cout << " C2 (Abs. Res) = " << C2_Value << "\n";
        //std::cout << " C3 (Rel. Res) = " << C3_Value << "\n";
        //std::cout << " C4 (Guar. Error) = " << C4_Value << "\n";

        if (Is_C1_Met || Is_C2_Met || Is_C3_Met) {
            std::cout << "Notification at Iteration " << k + 1 << ": Criterion ";

            bool first = true;
            if (Is_C1_Met) {
                std::cout << " (Relative Difference)";
                first = false;
            }
            if (Is_C2_Met) {
                if (!first) std::cout << " and ";
                std::cout << " (Absolute Residual)";
                first = false;
            }
            if (Is_C3_Met) {
                if (!first) std::cout << " and ";
                std::cout << " (Relative Residual)";
            }
            /*if (Is_C4_Met) {
                if (!first) std::cout << " and ";
                std::cout << "C4 (Guaranteed Error)";
            }*/
            std::cout << " has been met.\n";
        }

        if (Is_C1_Met && Is_C2_Met && Is_C3_Met) {
            std::cout << "All three criteria for convergence have been met after " << k + 1 << " iterations.\n";
            return x_next;
        }

        x_curr = x_next;
    }

    std::cout << "Maximum iterations (" << max_iter << ") reached. Did not converge within epsilon.\n";
    return x_curr;
}

template<typename T>
matrix<T> jacobi_method_elementwise(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T q,
    int max_iter, T epsilon)
{
    size_t n = A.rows;
    if (A.rows != A.cols || n != b.rows || b.cols != 1) {
        throw std::runtime_error("Invalid matrix/vector dimensions.");
    }

    matrix<T> x_curr = x0;
    matrix<T> x_next(n, 1, 0.0);

    std::cout << "Starting Element-wise Jacobi Method\n";
    const T epsilon_0 = static_cast<T>(1e-16);
    T C4_Threshold = ((static_cast<T>(1.0) - q) / q) * epsilon;
    matrix<T> Ax0_minus_b = A * x0 - b;
    T Norm_R0 = norminf(Ax0_minus_b);
    for (int k = 0; k < max_iter; ++k)
    {
        bool is_stable = true;

        for (size_t i = 0; i < n; ++i)
        {

            T a_ii = A.data[i * n + i];
            if (std::abs(a_ii) < 1e-12) {
                throw std::runtime_error("Zero diagonal element encountered.");
            }

            T sum = 0;

            for (size_t j = 0; j < n; ++j)
            {
                if (i != j) {
                    T a_ij = A.data[i * n + j];
                    T x_j_k = x_curr.data[j];
                    sum += a_ij * x_j_k;
                }
            }

            T b_i = b.data[i];


            x_next.data[i] = (b_i - sum) / a_ii;
        }
        matrix<T> residual_diff = x_next - x_curr;
        T diff_norm = norminf(residual_diff);

        matrix<T> R_next_vector = A * x_next - b;
        T R_next_norm = norminf(R_next_vector);

        T x_curr_norm = norminf(x_curr);


        T C1_Value = diff_norm / (x_curr_norm + epsilon_0);
        bool Is_C1_Met = (C1_Value <= epsilon);


        T C2_Value = R_next_norm;
        bool Is_C2_Met = (C2_Value <= epsilon);


        T C3_Value = (Norm_R0 < 1e-16) ? 0.0 : R_next_norm / Norm_R0;
        bool Is_C3_Met = (C3_Value <= epsilon);

        /*T C4_Value = diff_norm;
        bool Is_C4_Met = (C4_Value <= C4_Threshold);*/

        std::cout << "Iteration " << k + 1 << ":\n";
        std::cout << " C1 (Rel. Diff) = " << C1_Value << "\n";
        std::cout << " C2 (Abs. Res) = " << C2_Value << "\n";
        std::cout << " C3 (Rel. Res) = " << C3_Value << "\n";
        if (Is_C1_Met || Is_C2_Met || Is_C3_Met) {
            std::cout << "Notification at Iteration " << k + 1 << ": Criterion ";

            bool first = true;
            if (Is_C1_Met) {
                std::cout << " (Relative Difference)";
                first = false;
            }
            if (Is_C2_Met) {
                if (!first) std::cout << " and ";
                std::cout << " (Absolute Residual)";
                first = false;
            }
            if (Is_C3_Met) {
                if (!first) std::cout << " and ";
                std::cout << " (Relative Residual)";
            }
            /*if (Is_C4_Met) {
                if (!first) std::cout << " and ";
                std::cout << "C4 (Guaranteed Error)";
            }
            std::cout << " has been met.\n";*/
        }

        if (Is_C1_Met && Is_C2_Met && Is_C3_Met) {
            std::cout << "All three criteria for convergence have been met after " << k + 1 << " iterations.\n";
           // return x_next;
        }
        /*T max_diff = 0.0;
        for (size_t i = 0; i < n; ++i) {
            T diff = std::abs(x_next.data[i] - x_curr.data[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        std::cout << "Iteration " << k + 1 << ": ||x^(k+1) - x^k||_inf = " << max_diff << "\n";

        if (max_diff < epsilon) {
            std::cout << "Converged after " << k + 1 << " iterations.\n";
            return x_next;
        }*/


        x_curr = x_next;
    }

    std::cout << "Maximum iterations (" << max_iter << ") reached. Did not converge.\n";
    return x_curr;
}
template<typename T>
matrix<T> sor_method(const matrix<T>& A, const matrix<T>& b, const matrix<T>& x0, T omega, T norm_C, T norm_G2, int max_iter, T epsilon)
{
    size_t n = A.rows;

    matrix<T> x_curr = x0;
    matrix<T> x_prev(n, 1, 0.0);

    matrix<T> Ax0_minus_b = A * x0 - b;
    T Norm_R0 = norminf(Ax0_minus_b);

    const T epsilon_0 = static_cast<T>(1e-16);

    T C4_Threshold = ((static_cast<T>(1.0) - norm_C) / norm_G2) * epsilon;


    std::cout << "Starting SOR Method (omega = " << omega << ", ||C||=" << norm_C << ")\n";

    for (int k = 0; k < max_iter; ++k)
    {
        x_prev = x_curr;

        for (size_t i = 0; i < n; ++i)
        {
            T a_ii = A.data[i * n + i];
            if (std::abs(a_ii) < 1e-12) {
                throw std::runtime_error("Zero diagonal element encountered. Method fails.");
            }

            T sum_terms = 0;
            for (size_t j = 0; j < n; ++j)
            {
                if (i != j) {
                    T a_ij = A.data[i * n + j];
                    if (j < i) {
                        sum_terms += a_ij * x_curr.data[j];
                    }
                    else {
                        sum_terms += a_ij * x_prev.data[j];
                    }
                }
            }

            T b_i = b.data[i];
            T x_i_gs_approx = (b_i - sum_terms) / a_ii;
            x_curr.data[i] = (1.0 - omega) * x_prev.data[i] + omega * x_i_gs_approx;
        }

        matrix<T> diff_vector = x_curr - x_prev;
        T diff_norm = norminf(diff_vector);

        matrix<T> R_next_vector = A * x_curr - b;
        T R_next_norm = norminf(R_next_vector);

        T x_curr_norm = norminf(x_curr);

        T C1_Value = diff_norm / (x_curr_norm + epsilon_0);
        bool Is_C1_Met = (C1_Value <= epsilon);

        T C2_Value = R_next_norm;
        bool Is_C2_Met = (C2_Value <= epsilon);

        T C3_Value = (Norm_R0 < 1e-16) ? 0.0 : R_next_norm / Norm_R0;
        bool Is_C3_Met = (C3_Value <= epsilon);


       /* T C4_Value = diff_norm;
        bool Is_C4_Met = (C4_Value <= C4_Threshold);*/



        std::cout << "Iteration " << k + 1 << ":\n";
        std::cout << " C1 (Rel. Diff) = " << C1_Value << "\n";
        std::cout << " C2 (Abs. Res) = " << C2_Value << "\n";
        std::cout << " C3 (Rel. Res) = " << C3_Value << "\n";
       /* std::cout << " C4 (Guar. Error) = " << C4_Value  << "\n";*/

        if (Is_C1_Met || Is_C2_Met || Is_C3_Met) {
            std::cout << "Notification at Iteration " << k + 1 << ": Criterion ";

            bool first = true;
            if (Is_C1_Met) {
                std::cout << "C1 (Relative Difference)";
                first = false;
            }
            if (Is_C2_Met) {
                if (!first) std::cout << " and ";
                std::cout << "C2 (Absolute Residual)";
                first = false;
            }
            if (Is_C3_Met) {
                if (!first) std::cout << " and ";
                std::cout << "C3 (Relative Residual)";
                first = false;
            }
            
            std::cout << " has been met.\n";
        }

        if (Is_C1_Met && Is_C2_Met && Is_C3_Met) {
            std::cout << "All four criteria for convergence have been met after " << k + 1 << " iterations.\n";
            //return x_curr;
        }

    }

    std::cout << "Maximum iterations (" << max_iter << ") reached. Did not converge fully.\n";
    return x_curr;
}