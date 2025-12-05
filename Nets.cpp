#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>
#include "Nets.h"
#include "Matrices_template.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

template<typename T>
matrix<T> generate_uniform_grid_matrix(T a, T b, int n) {
    if (n < 1) {
        throw std::invalid_argument("n must be greater than 0.");
    }


    matrix<T> grid(1, n + 1);

    T h = (b - a) / n;

    for (int i = 0; i <= n; ++i) {
        T x_i = static_cast<T>(a + i * h);
        grid[i] = x_i;
    }
    return grid;
}

template<typename T>
matrix<T> generate_chebyshev_grid_matrix(T a, T b, int n) {
    if (n < 1) {
        throw std::invalid_argument("n must be greater than 0.");
    }

    // Создаем матрицу 1 x (n + 1) для хранения узлов
    matrix<T> grid(1, n + 1);

    T center = (a + b) / 2.0;
    T radius = (b - a) / 2.0;
    T factor = M_PI / (2.0 * (n + 1));

    for (int i = 0; i <= n; ++i) {
        // Аргумент косинуса: (2i + 1) * pi / (2(n + 1))
        T theta = (2.0 * i + 1.0) * factor;

        // Применяем формулу и приводим к типу T
        T x_i = static_cast<T>(center + radius * std::cos(theta));

        // Запись в матрицу
        grid[n - i] = x_i;
    }
    return grid;
}