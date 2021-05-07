#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <array>
#include <fstream>

#define METHODS 6
#define ITER_COUNT 30

template <typename T>
T function(T x) {
    return sin(x);
}

template <typename T>
T real_derivative(T x) {
    return cos(x);
}
template <typename T>
T backward_difference_2(T x, T h) {
    return (function(x) - function(x - h)) / h;
}

template <typename T>
T forward_difference_2(T x, T h) {
    return (function(x + h) - function(x)) / h;
}

template <typename T>
T central_difference_2(T x, T h) {
    return (function(x + h) - function(x - h)) / ((T)(2.0) * h);
}

template <typename T>
T forward_difference_3(T x, T h) {
    return (((T)(-3.0 / 2.0) * function(x)) + ((T)(2.0) * function(x + h)) - ((T)(1.0 / 2.0) * function(x + (T)2.0 * h))) / h;
}

template <typename T>
T backward_difference_3(T x, T h) {
    return (((T)(1.0 / 2.0) * function(x - (T)2.0 * h)) - ((T)(2.0) * function(x - h)) + ((T)(3.0 / 2.0) * function(x))) / h;
}

template <typename T>
void create_file(std::array<std::array<T, METHODS>, ITER_COUNT> errors, const std::string &file_name) {
    std::ofstream output(file_name);

    if (!output.good()){
        std::cout << "Blad podczas otwierania pliku" << std::endl;
        exit(1);
    }

    for (int i = 0; i < ITER_COUNT; i++) {
        for (int j = 0; j < METHODS; j++) {
            output << std::scientific << log10(errors[i][j]) << " ";
        }
        output << std::endl;
    }
    output.close();
}

template <typename T>
void calculate_errors(const std::string &file_name) {

    std::array<std::array<T, METHODS>, ITER_COUNT> errors;
    T h = 0.1;
    T start = 0;
    T end = M_PI / 2.0;
    T center = (start + end) / 2.0;

    for (int i = 0; i < ITER_COUNT; i++) {

        errors[i][0] = h;
        errors[i][1] = fabs(forward_difference_2(start, h) - real_derivative(start));
        errors[i][2] = fabs(forward_difference_3(start, h) - real_derivative(start));
        errors[i][3] = fabs(central_difference_2(center, h) - real_derivative(center));
        errors[i][4] = fabs(backward_difference_2(end, h) - real_derivative(end));
        errors[i][5] = fabs(backward_difference_3(end, h) - real_derivative(end));
        h *= (T)0.2;
        for (int j = 0; j < METHODS; j++) {
            printf("%.6e ", errors[i][j]);
        }
        std::cout << std::endl;
    }

    for (int i = 1; i < METHODS; i++) {
        std::cout << "Rzad dokladnosci:" << (log10(errors[3][i]) - log10(errors[2][i]))
                                                        / (log10(errors[3][0]) - log10(errors[2][0])) << std::endl;
    }
    std::cout << std::endl << std::endl;
    create_file(errors, file_name);
}

int main() {
    calculate_errors<float>("float.txt");
    calculate_errors<double>("double.txt");
    return 0;
}