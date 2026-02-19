#pragma once
//#include "Nets.cpp"
#include "Matrices_template.h"
//template<typename T>
//class Lagrange_polynom
//{
//public:
//    // Узлы интерполяции (x_i), 1xN. Должны быть доступны для конструктора
//    matrix<T> points;
//
//    // Коэффициенты полинома (разделенные разности f[x0], f[x0, x1], ...), 1xN
//    matrix<T> coefficients;
//
//    size_t N;               // Количество узлов
//
//    /**
//     * @brief Конструктор полинома Ньютона.
//     * * Вычисляет и сохраняет разделенные разности (коэффициенты).
//     * @param grid Матрица 1xN или Nx1, содержащая узлы интерполяции (x_i).
//     * @param f Указатель на функцию T(T), которую нужно аппроксимировать.
//     */
//    Lagrange_polynom(const matrix<T>& grid, T(*f)(T));
//
//    /**
//     * @brief Вычисляет значение полинома в точке x_eval по схеме Горнера.
//     * * @param x_eval Точка, в которой нужно вычислить полином.
//     * @return Значение полинома P(x_eval).
//     */
//    T operator()(T x_eval) const;
//
//    // --- Дополнительные объявления, если нужно, например, для отладки ---
//    // Объявляем, но не определяем методы доступа, если они нужны
//    const matrix<T>& get_coefficients() const { return coefficients; }
//};
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________

//template<typename T>
//class Lagrange_polynom
//{
//public:
//	matrix<T> points;
//	matrix<T> divided_differences;
//	size_t degree;
//	Lagrange_polynom() : points(), divided_differences(), degree(0) {}
//	Lagrange_polynom(matrix<T>& grid, T(*f)(T));
//	Lagrange_polynom(matrix<T>& grid, matrix<T>& values);
//
//	T operator()(T x)&;
//
//
//	~Lagrange_polynom();
//};
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________
//____________________________________Никитино__________________________________________


template<typename T>
class Lagrange_polynom
{
public:
    // Узлы интерполяции (x_i), 1xN.
    matrix<T> points;

    // Значения функции в узлах (y_i = f(x_i)), 1xN.
    matrix<T> values;

    size_t N;               // Количество узлов

    /**
     * @brief Конструктор полинома Лагранжа (стандартная форма).
     * * Сохраняет узлы и генерирует значения функции в узлах.
     * * Сложность: O(N).
     * @param grid Матрица узлов (x_i).
     * @param f Указатель на функцию T(T), которую нужно аппроксимировать.
     */
    Lagrange_polynom(const matrix<T>& grid, T(*f)(T));

    /**
     * @brief Вычисляет значение полинома в точке x_eval по стандартной формуле Лагранжа.
     * * Сложность: O(N^2).
     * @param x_eval Точка, в которой нужно вычислить полином.
     * @return Значение полинома P(x_eval).
     */
    T operator()(T x_eval) const;

    // Опционально: метод доступа к узлам
    const matrix<T>& get_points() const { return points; }
    const matrix<T>& get_values() const { return values; }
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