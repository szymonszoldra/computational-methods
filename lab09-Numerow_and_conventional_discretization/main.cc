#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

#define p 1.0
#define q 0.0
#define r -4.0
#define ALPHA 0.0
#define BETA 1.0
#define GAMMA -1.0
#define fi 0.0
#define psi 1.0
#define theta 0.0
#define interval_beginning 0.0
#define interval_end 1.0

#define EPSILON 1e-4

#define NUMEROW_PATH "numerow.txt"
#define CONV_PATH "conv.txt"
#define ERRORS_PATH "errors.txt"

bool save_numerow_to_file = true;
bool save_conv_to_file = true;

double analytical(double x) {
    return (exp(2.0 - 2.0 * x) - 4.0 * exp(4.0 - 2.0 * x) + 4.0 * exp(2.0 * x) - exp(2.0 + 2.0 * x) - x + x * exp(4.0)) / (4.0 - 4.0 * exp(4.0));
}

void thomas_algorithm(double* Lower, double* Diagonal, double* Upper, double* vector_b, double* vector_x, const int steps) {
    for (int i = 1; i < steps; i++) {
        Diagonal[i] = Diagonal[i] - (Upper[i] * Lower[i - 1]) / Diagonal[i - 1];
    }

    for (int i = 1; i < steps; i++) {
        vector_b[i] = vector_b[i] - (Upper[i] * vector_b[i - 1]) / Diagonal[i - 1];
    }

    for (int i = steps - 1; i >= 0; i--) {
        if (i == steps -1) {
            vector_x[i] = vector_b[i] / Diagonal[i];
            continue;
        }
        vector_x[i] = (vector_b[i] - Lower[i] * vector_x[i + 1]) / Diagonal[i];
    }
}

double max_error(double* errors, const int steps) {
    double current_max = fabs(errors[0]);

    for (int i = 1; i < steps; i++) {
        current_max = std::max(current_max, fabs(errors[i]));
    }

    delete[] errors;
    return current_max;
}

double discretization(const double h, const int steps, const bool numerow) {
    double *Lower     = new double[steps];
    double *Diagonal  = new double[steps];
    double *Upper     = new double[steps];
    double *vector_b  = new double[steps];
    double *vector_x  = new double[steps];
    double *errors    = new double[steps];

    Lower[0] = 0.0;
    Diagonal[0] = BETA - ALPHA / h;
    Upper[0] = ALPHA / h;
    vector_b[0] = -GAMMA;

    if (numerow) {
        for (int i = 1; i < steps - 1; i++) {
            Lower[i] = p / (h * h) + r / 12.0;
            Diagonal[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
            Upper[i] = p / (h * h) + r / 12.0;
            vector_b[i] = (i * h - h) / 12.0 - (10.0 / 12.0) * -1.0 * (i * h) + (i * h + h) / 12.0;
        }
    } else {
        for (int i = 1; i < steps - 1; i++) {
            Lower[i] = p / (h * h) - q / (2.0 * h);
            Diagonal[i] = (-2.0 * p) / (h * h) + r;
            Upper[i] = p / (h * h) - q / (2.0 * h);
            vector_b[i] = i * h;
        }
    }

    Lower[steps - 1] = -fi / h;
    Diagonal[steps - 1] = -fi / h + psi;
    vector_b[steps - 1] = -theta;
    Upper[steps - 1] = 0.0;

    thomas_algorithm(Lower, Diagonal, Upper, vector_b, vector_x, steps);

    // I only needed to write to file once when there were 100 steps for each conv and Numerow
    if (save_conv_to_file || save_numerow_to_file) {
        if (numerow) {
            save_numerow_to_file = false;
        } else {
            save_conv_to_file = false;
        }

        std::fstream output;
        output.open(numerow ? NUMEROW_PATH : CONV_PATH, std::fstream::out);
        double file_steps = interval_beginning;
        for (int i = 0; i < steps; i++) {
            output << file_steps << " " << vector_x[i] << " " << analytical(file_steps) << std::endl;
            file_steps += h;
        }
    }

    double xi = interval_beginning;
    for (int i = 0; i < steps; i++) {
        errors[i] = fabs(vector_x[i] - analytical(xi));
        xi += h;
    }

    delete[] Lower;
    delete[] Diagonal;
    delete[] Upper;
    delete[] vector_b;
    delete[] vector_x;

    return max_error(errors, steps);
}

int main() {
    int steps = 100;
    double error_conv = 0;
    double error_numerow = 0;
    double distance = interval_end - interval_beginning;
    double h = distance / (steps - 1);

    std::fstream errors;
    errors.open(ERRORS_PATH, std::fstream::out);

    while (h > EPSILON) {
        h = distance / (steps - 1);

        // Both conventional and Numerow discretization only differ slightly in one for loop
        // so there is no need for them to be two separate functions. 3rd argument named "numerow"
        // takes care for that.
        error_conv = log10(discretization(h, steps, false));
        error_numerow = log10(discretization(h, steps, true));
        errors << log10(h) << " " << error_conv << " " << error_numerow << std::endl;
        steps += 100;
    }

    errors.close();
    return 0;
}
