#pragma once
#include <iostream>
#include "Matrices_template.cpp"
#include "Nelinear.h"

template<typename T>
matrix<T> localize_roots(T a, T b, T step, T(*f)(T)) {
    // Используем matrix<double> как динамический массив, в котором храним
    // пары значений [x_left, x_right]
    matrix<T> localization_intervals(0, 0); // Инициализация пустым массивом

    if (step <= 0 || a >= b) {
        std::cerr << "Ошибка: Шаг должен быть > 0, и a < b." << std::endl;
        return localization_intervals;
    }

    size_t current_size = 0;
    T x_i = a;
    T f_i = f(x_i);

    while (x_i < b) {
        T x_i_plus_1 = x_i + step;
        // Убедимся, что не выходим за b
        if (x_i_plus_1 > b) {
            x_i_plus_1 = b;
        }

        // Если шаг слишком мал (например, x_i = b), выходим
        if (x_i_plus_1 <= x_i) {
            break;
        }

        T f_i_plus_1 = f(x_i_plus_1);

        // Проверка условия смены знака: f(x_i) * f(x_{i+1}) < 0
        if (f_i * f_i_plus_1 < 0) {
            // Найден интервал локализации: [x_i, x_i_plus_1]

            // 1. Создаем новый массив с увеличенным размером (в 2 раза больше для пары [x_left, x_right])
            size_t new_size = current_size + 2;
            matrix<T> temp_matrix(1, new_size);

            // 2. Копируем старые данные (если есть)
            for (size_t k = 0; k < current_size; ++k) {
                temp_matrix[k] = localization_intervals[k];
            }

            // 3. Добавляем новый интервал
            temp_matrix[current_size] = x_i;
            temp_matrix[current_size + 1] = x_i_plus_1;

            // 4. Обновляем localization_intervals (используем swap для эффективности)
            std::swap(localization_intervals.data, temp_matrix.data);
            localization_intervals.rows = 1;
            localization_intervals.cols = new_size;
            /*localization_intervals.size_t current_size = new_size;*/

            // Важно: в реальном коде необходимо корректно использовать swap
            // или operator= для matrix<T> во избежание утечек памяти
            // (здесь используется упрощенный подход)

            current_size = new_size;
            localization_intervals.rows = 1;
            localization_intervals.cols = current_size;

        }
        else if (f_i == 0.0) {
            // Обнаружен корень на границе интервала, добавляем [x_i, x_i]
            // (Повторяем процесс расширения массива, как выше)
            size_t new_size = current_size + 2;
            matrix<T> temp_matrix(1, new_size);
            for (size_t k = 0; k < current_size; ++k) {
                temp_matrix[k] = localization_intervals[k];
            }
            temp_matrix[current_size] = x_i;
            temp_matrix[current_size + 1] = x_i;

            std::swap(localization_intervals.data, temp_matrix.data);
            current_size = new_size;
            localization_intervals.rows = 1;
            localization_intervals.cols = current_size;
        }

        // Переходим к следующему интервалу
        x_i = x_i_plus_1;
        f_i = f_i_plus_1; // Переиспользуем уже вычисленное значение
    }

    return localization_intervals;
}

template<typename T>
T bisection_method(T a, T b,T(*f)(T),T tolerance, int& iterations) {
    iterations = 0;
    T c; // Середина интервала

    // Проверка, что на концах интервала функция имеет разные знаки
    if (f(a) * f(b) >= 0) {
        std::cerr << "Error: The bisection method is not applicable. The function does not change sign over the interval." << std::endl;
        return std::nan(""); // Возвращаем "не число"
    }

    while ((b - a) >= tolerance) {
        iterations++;
        c = (a + b) / 2.0;

        if (f(c) == 0.0) {
            return c; // Точное совпадение
        }
        else if (f(c) * f(a) < 0) {
            b = c; // Корень в левой половине
        }
        else {
            a = c; // Корень в правой половине
        }
    }

    return (a + b) / 2.0;
}

template<typename T>
T newton_method(T a, T b,T x0, T(*f)(T), T(*derivative_f)(T), T tolerance,  int max_iterations, int& iterations) {
    iterations = 0;
    T x_curr = x0;
    T x_next;
    T f_x;
    T f_prime_x;
    for (iterations = 0; iterations < max_iterations; ++iterations) {
        f_x = f(x_curr);
        f_prime_x = derivative_f(x_curr);

        // Проверка на деление на ноль (производная близка к нулю)
        if (std::abs(f_prime_x) < 1e-10) {
            std::cerr << "Error: The derivative is close to zero. Newton's method does not converge." << std::endl;
            return std::nan("");
        }

        // Формула метода Ньютона
        T delta_x = x_curr - f_x / f_prime_x;
        T alpha = 1.0;
        T x_next;
        const T min_alpha = 1e-12; // Минимально допустимый шаг
        
        // --- МЕХАНИЗМ LINE SEARCH С ПРОВЕРКОЙ ОО ---
        while (true) {
            x_next = x_curr + alpha * delta_x;
            
            // 1. Проверка области определения
            if (a<=x_next<=b) {
                // Точка внутри ОО. Шаг принят.
                break;
            }
            
            // 2. Уменьшение шага
            alpha /= 2.0;
            
            // 3. Проверка на минимальный шаг (если alpha слишком маленькое)
            if (alpha < min_alpha) {
                std::cerr << "Error: Cannot find a step to stay within the domain." << std::endl;
                return x_curr; // Возвращаем последнее рабочее приближение
            }
        }
        // Критерий остановки: |x_{k+1} - x_k| < tolerance
        if (std::abs(x_next - x_curr) < tolerance) {
            return x_next;
        }

        x_curr = x_next;
    }

    std::cerr << "Warning: Newton's method did not work in" << max_iterations << " iterations." << std::endl;
    return x_curr; // Возвращаем последнее приближение
}

template<typename T>
void find_roots(T a, T b, T(*f)(T), T(*derivative_f)(T), T tolerance){
    int bisection_iters = 0;
    int newton_iters = 0;
    const int MAX_NEWTON_ITERS = 100; // Ограничение для Ньютона
    T x0_newton = (a + b) / 2.0; // Начальное приближение для Ньютона - середина интервала

    std::cout << "\n--- Interval: [" << std::fixed << a << ", " << b << "] ---" << std::endl;

    // --- Метод Бисекции ---
    T root_bisection = bisection_method(a, b, f,tolerance, bisection_iters);
    if (!std::isnan(root_bisection)) {
        std::cout << "Bisection Method :\n";
        std::cout << "  Root: " << std::fixed << root_bisection << "\n";
        std::cout << "  Iterations: " << bisection_iters << "\n";
    }

    // --- Метод Ньютона ---
    T root_newton = newton_method(a,b,x0_newton,f, derivative_f, tolerance, MAX_NEWTON_ITERS, newton_iters);
    if (!std::isnan(root_newton)) {
        std::cout << "Newton's method:\n";
        std::cout << "  Root: " << std::fixed << root_newton << "\n";
        std::cout << "  Iterations: " << newton_iters << "\n";
    }
}

template<typename T>
// Численно вычисляет матрицу Якоби 2x2.
matrix<T> numerical_jacobian(const matrix<T>& x_curr, matrix<T>(*f_system)(const matrix<T>&) ) {
    T H = 1e-10;
    
    matrix<T> J(2, 2);
    matrix<T> f_curr = f_system(x_curr); // Вектор f(x_curr)

    // Копирование для работы с x1 и x2
    matrix<T> x_temp = x_curr;

    // --- Первый столбец: Производные по x1 ---

    // x1 + h
    x_temp(0, 0) = x_curr(0, 0) + H;
    matrix<T> f_plus_h1 = f_system(x_temp);

    // J11 = (f1(x1+h, x2) - f1(x1, x2)) / h
    J(0, 0) = (f_plus_h1(0, 0) - f_curr(0, 0)) / H;
    // J21 = (f2(x1+h, x2) - f2(x1, x2)) / h
    J(1, 0) = (f_plus_h1(1, 0) - f_curr(1, 0)) / H;

    // Сброс x_temp(0, 0) до исходного x1
    x_temp(0, 0) = x_curr(0, 0);

    // --- Второй столбец: Производные по x2 ---

    // x2 + h
    x_temp(1, 0) = x_curr(1, 0) + H;
    matrix<T> f_plus_h2 = f_system(x_temp);

    // J12 = (f1(x1, x2+h) - f1(x1, x2)) / h
    J(0, 1) = (f_plus_h2(0, 0) - f_curr(0, 0)) / H;
    // J22 = (f2(x1, x2+h) - f2(x1, x2)) / H
    J(1, 1) = (f_plus_h2(1, 0) - f_curr(1, 0)) / H;

    return J;
}

template<typename T>
matrix<T> solve_newton_step(const matrix<T>& f_vec, const matrix<T>& J) {
    T DET_TOL = 0.0000001;

    T R1 = -f_vec(0, 0);
    T R2 = -f_vec(1, 0);

    // Элементы матрицы Якоби
    T a = J(0, 0), b = J(0, 1), c = J(1, 0), d = J(1, 1);
    T detJ = a * d - b * c;

    if (std::abs(detJ) < DET_TOL) {
        std::cerr << "Error: Jacobian matrix is singular (det=0)." << std::endl;
        // Возвращаем вектор с NaN, как индикатор ошибки
        return matrix<T>(2, 1, std::numeric_limits<T>::quiet_NaN());
    }

    // Решение по Правилу Крамера
    T dx1 = (R1 * d - b * R2) / detJ;
    T dx2 = (a * R2 - R1 * c) / detJ;

    matrix<T> delta_x(2, 1);
    delta_x(0, 0) = dx1;
    delta_x(1, 0) = dx2;

    return delta_x;
}
template<typename T>
T norm_L2(const  matrix<T>& v) {
    return std::sqrt(v(0, 0) * v(0, 0) + v(1, 0) * v(1, 0));
}

template<typename T>
matrix<T> newton_system_2x2_matrix(
    T a1,
    T b1,
    T a2,
    T b2,
    const matrix<T>& x0,
    matrix<T>(*f_system)(const matrix<T>&),
    /*DomainPredicate is_defined,*/
    T tolerance,
    int max_iterations,
    int& iterations)
{
    iterations = 0;
    matrix<T> x_curr = x0;

    if (!(a1<=x_curr(0,0)<=b1&& a2 <= x_curr(1, 0) <= b2)) {
        std::cerr << "Error: Initial guess x0 is outside the domain." << std::endl;

        iterations = 100;
        return matrix<T>(2, 1, std::numeric_limits<T>::quiet_NaN());
    }

    const T min_alpha = 1e-12; // Минимально допустимый шаг

    for (iterations = 0; iterations < max_iterations; ++iterations) {

        matrix<T> f_x = f_system(x_curr);

        // 1. Критерий остановки по норме вектора F
        /*if (norm_L2(f_x) < tolerance) {
            iterations = max_iterations + 3;
            return x_curr;
        }*/

        // 2. Вычисление Якоби и шага
        matrix<T> J = numerical_jacobian(x_curr, f_system);
        matrix<T> delta_x = solve_newton_step(f_x, J);

        // Проверка на ошибку сингулярности
        if (std::isnan(delta_x(0, 0))) {
            iterations = 100;
            return x_curr;
        }

        // 3. МЕХАНИЗМ LINE SEARCH С ПРОВЕРКОЙ ОО
        T alpha = 1.0;
        matrix<T> x_next;

        while (true) {
            // x_next = x_curr + alpha * delta_x
            x_next = x_curr + delta_x * alpha;

            if (a1 <= x_curr(0, 0) <= b1 && a2 <= x_next(0, 1) <= b2) {
                break; // Точка внутри ОО. Шаг принят.
            }

            // Уменьшение шага
            alpha /= 2.0;

            if (alpha < min_alpha) {
                std::cerr << "Warning: Cannot find a step to stay within the domain. Stalled." << std::endl;
                iterations = max_iterations + 2;
                return x_curr;
            }
        }
        // КОНЕЦ LINE SEARCH

        // 4. Критерий остановки по изменению X
        // norm(x_next - x_curr) < tolerance
        if (norm_L2(x_next - x_curr) < tolerance) {

            return x_next;
        }

        x_curr = x_next;
    }

    std::cerr << "Warning: Maximum iterations reached." << std::endl;
    iterations = max_iterations + 1;
    return x_curr;
}
template<typename T>
int run_newton_for_map(T a1,T b1,T a2, T b2,
    const matrix<T>& x0,
    matrix<T>(*f_system)(const matrix<T>&),
    T tolerance,
    int max_iterations)
{
    if (!(a1 <= x0(0, 0) <= b1 && a2 <= x0(1, 0) <= b2)) {
        return max_iterations + 2;
    }

    matrix<T> x_curr = x0;
    T min_alpha = 1e-12;

    for (int iterations = 0; iterations < max_iterations; ++iterations) {

        matrix<T> f_x = f_system(x_curr);

        // Критерий остановки: f(x) близка к 0
        if (norm_L2(f_x) < tolerance) {
            return iterations;
        }

        // Вычисление Якоби и шага
        matrix<T> J = numerical_jacobian(x_curr, f_system);
        matrix<T> delta_x = solve_newton_step(f_x, J);

        // Проверка на ошибку сингулярности
        if (std::isnan(delta_x(0, 0))) {
            return max_iterations+3;
        }

        // LINE SEARCH И ПРОВЕРКА ОО
        T alpha = 1.0;
        matrix<T> x_next;

        while (true) {
            // x_next = x_curr + delta_x * alpha; (Использует operator+ и operator*(T k))
            x_next = x_curr + delta_x * alpha;

            if (!(a1 <= x_next(0, 0) <= b1 && a2 <= x_next(1, 0) <= b2)) {
                break; // Точка внутри ОО. Шаг принят.
            }
            alpha /= 2.0;
            if (alpha < min_alpha) {
                return max_iterations+1; // Невозможно сделать шаг, оставаясь в ОО
            }
        }

        // Критерий остановки 2: Изменение x_i мало
        if (norm_L2(x_next - x_curr) < tolerance) {
            return iterations;
        }

        x_curr = x_next;
    }

    return max_iterations + 1;
}