#include <iostream>
#include <cmath>
#include <array>

#define MAX_ITER 100
#define TOLF 1e-8
#define TOLX 1e-8

void print_vec(const std::array<double, 3> &a) {
    std::cout << "[ ";
    for (const double i: a) {
        printf("%.16f ", i);
    }
    std::cout << "]" << std::endl;
}

bool should_stop_iterate (const std::array<double, 3> &en, const std::array<double, 3> &residue) {
    int counter = 0;

    for (const double i: en) {
        if(i <= TOLX)
            counter++;
    }

    for (const double i: residue) {
        if(fabs(i) <= TOLF)
            counter++;
    }

    // The residue rule and the error estimator rule for each unknown must be satisfied
    if (counter == 6)
        return true;

    return false;
}

std::array<double, 3> correction(const std::array<double, 3> &a) {
    double x = a[0];
    double y = a[1];
    double z = a[2];

    // The correction for functions from line 50 (calculated manually)
    return std::array<double, 3>{
            (2.0 * x * x * y - y * y + x * x - 1.0) / (4.0 * x * y + 2.0 * x),
            (y * y + y - 1.0) / (2.0 * y + 1.0),
            (z * z - 1.0) / (2.0 * z)
    };
}

std::array<double, 3> functions(const std::array<double, 3> &a) {
    double x = a[0];
    double y = a[1];
    double z = a[2];

    return std::array<double, 3>{
            x * x + y * y + z * z - 2.0,
            x * x + y * y - 1.0,
            x * x - y
    };
}

int main() {
    std::array<double, 3> x_n = {1.0, 1.0, 1.0};
    print_vec(x_n);

    std::array<double, 3> x_n1;
    std::array<double, 3> en;
    std::array<double, 3> residue;
    std::array<double, 3> next_correction;


    for (int i = 1; i <= MAX_ITER; i++) {
        next_correction = correction(x_n);

        x_n1[0] = x_n[0] - next_correction[0];
        x_n1[1] = x_n[1] - next_correction[1];
        x_n1[2] = x_n[2] - next_correction[2];

        en[0] = fabs(x_n[0] - x_n1[0]);
        en[1] = fabs(x_n[1] - x_n1[1]);
        en[2] = fabs(x_n[2] - x_n1[2]);

        residue = functions(x_n1);


        x_n[0] = x_n1[0];
        x_n[1] = x_n1[1];
        x_n[2] = x_n1[2];

        printf("Obecna iteracja: %2d \n", i);
        printf("Zmienne: ");
        print_vec(x_n);
        printf("Estymatory: ");
        print_vec(en);
        printf("Residua: ");
        print_vec(residue);

        std::cout << std::endl;

        if(should_stop_iterate(en, residue))
            break;
    }

    return 0;
}
