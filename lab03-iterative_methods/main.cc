#include "functions.h"

#include <cmath>
#include <iostream>
#include <functional>

#define MAX_ITER 50
#define TOLF 1e-8
#define TOLX 1e-8

void picard_method(const std::function<double(double)> &function,
                   const std::function<double(double)> &function_picard,
                   const std::function<double(double)> &function_picard_derivative,
                   double x,
                   const char* name
) {

    std::cout << std::endl << "Metoda Picarda:"<< std::endl;

    if (fabs(function_picard_derivative(x)) >= 1){
        std::cout << "["<< name <<"] Wartość bezwzględna pochodnej Phi(x) nie jest mniejsza od 1, Rozbieżność"
                  << std::endl;
        return;
    }

    double en = x;
    double residue = x;
    double next_approx = x;

    for (int i = 1; i <= MAX_ITER; i++) {
        next_approx = function_picard(next_approx);
        en = fabs(next_approx - x);
        x = next_approx;
        residue = fabs(function(x));
        printf("[%s] %2d xn %.16f Estymator %.16f F(xn) %.16f\n", name, i, next_approx, en, residue);
        if ((residue <= TOLF) && (en <= TOLX)) {
            break;
        }
    }
}

void bisection_method(const std::function<double(double)> &function, double a, double b, const char* name) {
    std::cout << std::endl << "Metoda Bisekcji:" << std::endl;

    if (function(a) * function(b) > 0) {
        std::cout << "Podano zły przedział!" << std::endl;
        return;
    }

    double left = a;
    double right = b;
    double center = (left + right) / 2.0;
    double en = (right - left) / 2.0;
    double residue = fabs(function(center));

    for (int i = 1; i <= MAX_ITER; i++) {

        printf("[%s] %2d Left %.16f Right %.16f Center %.16f Estymator %.16f F(center) %.16f\n",
               name, i, left, right, center, en, residue);

        if (fabs(en) <= TOLX && fabs(residue) <= TOLF || function(center) == 0) {
            break;
        }

        if (function(left) * function(center) < 0) {
            right = center;
        } else {
            left = center;
        }

        center = (left + right) / 2.0;
        en = (right - left) / 2.0;
        residue = fabs(function(center));
    }
}

void newton_method(const std::function<double(double)> &function,
                   const std::function<double(double)> &function_derivative,
                   double x,
                   const char* name
) {
    std::cout << std::endl << "Metoda Newtona:" << std::endl;

    double x_n = x;
    double x_n1;
    double en;
    double residue = fabs(function(x_n));
    for (int i = 1; i <= MAX_ITER; i++) {
        x_n1 = x_n - function(x_n) / function_derivative(x_n);
        en = fabs(x_n - x_n1);
        x_n = x_n1;
        residue = fabs(function(x_n));

        printf("[%s] %2d x_n+1 %.16f Estymator %.16f F(x_n+1) %.16f\n", name, i, x_n, en, residue);
        if (residue <= TOLF && en <= TOLX) {
            break;
        }
    }
}

void secant_method(const std::function<double(double)> &function, double x_n, double x_n1, const char* name) {
    double x_n2 = x_n1 - function(x_n1) / ((function(x_n1) - function(x_n)) / (x_n1 - x_n));
    double en = fabs(x_n1 - x_n);
    double residue = fabs(function(x_n2));
    std::cout << std::endl << "Metoda Siecznych:" << std::endl;
    for (int i = 1; i <= MAX_ITER; i++) {
        printf("[%s] %2d x_n+2 %.16f Estymator %.16f F(x_n+2) %.16f\n",name, i, x_n2, en, residue);
        x_n = x_n1;
        x_n1 = x_n2;

        x_n2 = x_n1 - function(x_n1) / ((function(x_n1) - function(x_n)) / (x_n1 - x_n));
        if (residue <= TOLF && en <= TOLX) {
            break;
        }
        en = fabs(x_n1 - x_n);
        residue = fabs(function(x_n2));
    }
}

int main() {
    const char* name = "sin^2(x/4) - x = 0";
    picard_method(sin_function, sin_function_picard, sin_function_picard_derivative, 1.0, name);
    bisection_method(sin_function, -1.0, 1.0, name);
    newton_method(sin_function, sin_function_derivative, 1.0, name);
    secant_method(sin_function, 1.0, 0.5, name);

    name = "tan(2x)-x-1=0";
    picard_method(tan_function, tan_function_picard, tan_function_picard_derivative, 0.5, name);
    bisection_method(tan_function, 2.6, 3.9 , name);
    newton_method(tan_function, tan_function_derivative, 3.9, name);
    secant_method(tan_function, 3.9, 2.6, name);
    return 0;
}