#include <array>
#include <cmath>
#include <iostream>
#include <fstream>

// interval [a, b]
#define a -1.0
#define b 1.0

#define NODES 20
#define NUMBER_OF_STEPS 200
#define STEP (b - a) / NUMBER_OF_STEPS
#define RESULTS_PATH "results.txt"

double polynomial(const double x) {
    return 1.0 / (1.0 + 10.0 * x * x * x * x * x * x);
}

void Lagrange_Newton(std::array<std::array<double, 2>, NODES>& newton_array) {
    const int n = NODES - 1;
    int i = 0;
    for (double x = a; x <= b; x += ((b - a) / n), i++) {
        newton_array[i][0] = x;
        newton_array[i][1] = polynomial(x);
    }

    for (int j = 1; j < n + 1; j++) {
        for (i = n; i >= j; i--) {
            newton_array[i][1] = (newton_array[i][1] - newton_array[i - 1][1]) /
                                (newton_array[i][0] - newton_array[i - j][0]);
        }
    }
}

void Lagrange_Czebyszew(std::array<std::array<double, 2>, NODES>& czebyszew_array) {
    const int n = NODES - 1;

    for (int i = 0; i <= n; i++) {
        czebyszew_array[i][0] = cos(((2.0 * i + 1.0) * M_PI) / (2.0 * n + 2));
        czebyszew_array[i][1] = polynomial(czebyszew_array[i][0]);
    }

    for (int j = 1; j < n + 1; j++) {
        for (int i = n; i >= j; i--) {
            czebyszew_array[i][1] = (czebyszew_array[i][1] - czebyszew_array[i - 1][1]) /
                                    (czebyszew_array[i][0] - czebyszew_array[i - j][0]);
        }
    }
}

double calculate_values(const double x, const std::array<std::array<double, 2>, NODES>& arr) {
    const int n = NODES - 1;
    double value = arr[n][1];

    for (int i = n - 1; i >= 0; i--) {
        value *= (x - arr[i][0]);
        value += arr[i][1];
    }

    return value;
}

int main() {
    // 2d array with shape NODES x 2
    // on index 0 x values
    // on index 1 corresponding y values
    std::array<std::array<double, 2>, NODES> newton_array;
    std::array<std::array<double, 2>, NODES> czebyszew_array;

    Lagrange_Newton(newton_array);
    Lagrange_Czebyszew(czebyszew_array);

    std::fstream results;
    results.open(RESULTS_PATH, std::fstream::out);
    for (double x = a; x < b; x += STEP) {
        results << x
                << " "
                << calculate_values(x, czebyszew_array)
                << " "
                << calculate_values(x, newton_array)
                << std::endl;
    }
    results.close();
    return 0;
}