
#pragma once
#include <iostream>
#include "Matrices_template.h"

template<typename T>
matrix<T> localize_roots(T a, T b, T step, T(*f)(T));