#include <cmath>
#include <iostream>
#include <fstream>
#include <functional>

#define N 1000
#define EPSILON 1e-15

#define ERRORS_PATH "../errors.txt"
#define STABLE_PATH "../stable.txt"
#define UNSTABLE_PATH "../unstable.txt"

// Bezpośrednia metoda Eulera
double BME(const double h, const double t_m) {
    double y = 0.0;
    for (double t = 0.0; t < t_m; t += h) {
        y += h * (-((10.0 * t * t + 20.0) / (t * t + 1.0) * (y - 1.0)));
    }
    return y;
}

// Pośrednia metoda Eulera
double PME(const double h, const double t_m) {
    double y = 0.0;
    double temp;
    for (double t = 0.0; t < t_m; t += h) {
        temp = (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
        y = (y + h * temp) / (1.0 + h * temp);
    }
    return y;
}

// Metoda trapezów
double PMT(const double h, const double t_m) {
    double y = 0.0;
    double temp1;
    double temp2;
    for (double t = 0.0; t < t_m; t += h) {
        temp1 = (10.0 * t * t + 20.0) / (t * t + 1.0);
        temp2 = (10.0 * (t + h) * (t + h) + 20.0) / ((t + h) * (t + h) + 1.0);
        y = (y - (h * (temp1 * (y - 1.0) - temp2)) * 0.5) / (1.0 + ((h * temp2) * 0.5));
    }
    return y;
}

double analytical(double t) {
    return 1.0 - exp(-10.0 * (t + atan(t)));
}

double max_error(std::function<double(double, double)> method, const double h) {
    double t = h;
    double current_error = std::fabs(analytical(t) - method(h, t));

    for (int i = 0; i < N; i++) {
        current_error = std::max(std::fabs(analytical(t) - method(h, t)), current_error);
        t += h;
    }

    return log10(current_error);
}

int main() {
    double step = 0.01;

    std::ofstream errors;
    errors.open(ERRORS_PATH);
    for (double i = 20; i > EPSILON; i /= 1.5) {
        errors << log10(i);
        errors << " ";
        errors << max_error(BME, i);
        errors << " ";
        errors << max_error(PME, i);
        errors << " ";
        errors << max_error(PMT, i);
        errors << std::endl;
    }
    errors.close();

    std::ofstream stable;
    stable.open(STABLE_PATH);
    for (double t = 0.0; t < 5.0; t += step) {
        stable << t;
        stable << " ";
        stable << analytical(t);
        stable << " ";
        stable << BME(step, t);
        stable << " ";
        stable << PME(step, t);
        stable << " ";
        stable << PMT(step, t);
        stable << std::endl;
    }
    stable.close();

    double unstableStep = 0.2;
    std::ofstream unstable;
    unstable.open(UNSTABLE_PATH);
    for (double t = 0.0; t < 5.0; t += unstableStep) {
        unstable << t;
        unstable << " ";
        unstable << analytical(t);
        unstable << " ";
        unstable << BME(unstableStep, t);
        unstable << std::endl;
    }
    unstable.close();
    return 0;
}