#ifndef LAB3_FUNCTIONS_H
#define LAB3_FUNCTIONS_H
#include <cmath>
double sin_function(double x) {
    return sin(x / 4.0) * sin(x / 4.0) - x;
}

double sin_function_derivative(double x) {
    return 0.25 * sin(x / 2.0) - 1.0;
}

double sin_function_picard(double x) {
    return sin(x / 4.0) * sin(x / 4.0);
}

double sin_function_picard_derivative(double x) {
    return 0.25 * sin(x / 2.0);
}

double tan_function(double x) {
    return tan(2.0 * x) - x - 1.0;
}

double tan_function_derivative(double x) {
  return -1.0 + 2.0 / (cos(2.0 * x) * cos(2.0 * x));
}

double tan_function_picard(double x) {
    return tan(2.0 * x) - 1.0;
}

double tan_function_picard_derivative(double x) {
  return 2.0 / (cos(2.0 * x) * cos(2.0 * x));
}
#endif  // LAB3_FUNCTIONS_H