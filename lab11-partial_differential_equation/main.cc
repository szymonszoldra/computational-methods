#include <iostream>
#include <cmath>
#include <fstream>

#include "lib/calerf/calerf.h"

// consts

const int D = 1;
const double t_min = 0.0;
const double t_max = 2.0;

const double a = 6.0 * sqrt(D * t_max);

const double x_start = 0.0;
const double x_end = a;

const double KMB_LAMBDA = 0.4;
const double LAASONEN_LAMBDA = 1.0;
const double h = 0.1;

typedef double** matrix_t;
typedef double* vector_t;

struct Errors {
    matrix_t errors_matrix;
    vector_t max_error;
};

struct Steps {
    vector_t t_steps;
    vector_t x_steps;
};


namespace utils {
    void cleanup(matrix_t analytical, matrix_t computational, Errors err, Steps steps, const int n) {
            for (int i = 0; i < n; i++) {
                delete[] analytical[i];
                delete[] computational[i];
                delete[] err.errors_matrix[i];
            }
            delete[] analytical;
            delete[] computational;
            delete[] err.errors_matrix;
            delete[] err.max_error;
            delete[] steps.t_steps;
            delete[] steps.x_steps;
    }

    namespace vector {
        void write(vector_t vec, const int m, const char* name) {
            std::fstream file(name, std::ios::out);
            for (int i = 0; i < m; i++) {
                file << vec[i] << std::endl;
             }
            file.close();
        }

        double max_error(vector_t vec, const int m) {
            double current_max = fabs(vec[0]);
            for (int i = 1; i < m; i++) {
                current_max = std::max(current_max, fabs(vec[i]));
            }
            return current_max;
        }
    }

    namespace matrix {
        void write(matrix_t mtx, const int n, const int m, const char* name) {
            std::fstream file(name, std::ios::out);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) {
                    file << mtx[i][j] << " ";
                }
                file << std::endl;
            }
            file.close();
        }

        matrix_t get_matrix_with_conditions(const int n, const int m) {
            matrix_t matrix = new vector_t[n];
            for (int i = 0; i < n; i++) {
                matrix[i] = new double[m];
            }

            // Initial condition
            for (int i = 0; i < m; i++) {
                matrix[0][i] = 1.0;
            }

            // Boundary conditions
            for (int i = 0; i < n; i++) {
                matrix[i][0] = 0.0;
                matrix[i][m - 1] = 1.0;
            }

            return matrix;
        }
    }
}

double get_delta_t(const double lambda, const double h, const double D) {
    return (lambda * h * h) / D;
}

matrix_t get_analytic_solution(const int n, const int m, const double h, const double delta_t) {
    matrix_t matrix = new vector_t[n];
    for (int i = 0; i < n; i++) {
        matrix[i] = new double[m];
    }

    double x = x_start;
    double t = t_min;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = calerf::ERFL(x / (2 * sqrt(D * t)));
            x += h;
        }
        x = x_start;
        t += delta_t;
    }

    return matrix;
}

Errors get_error(matrix_t analytic, matrix_t computational_solution, const int n, const int m) {
    matrix_t errors = new vector_t[n];
    for (int i = 0; i < n; i++) {
        errors[i] = new double[m];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            errors[i][j] = fabs(analytic[i][j] - computational_solution[i][j]);
        }
    }

    vector_t max_error = new double[n];

    for (int i = 0; i < n; i++) {
        max_error[i] = utils::vector::max_error(errors[i], m);
    }

    return {
            errors, max_error
    };
}

Steps get_steps(double delta_t, int n, int m) {
    vector_t t_steps = new double[n];
    vector_t x_steps = new double[m];

    double x = x_start;
    double t = t_min;

    for (int i = 0; i < n; i++) {
        t_steps[i] = t;
        t += delta_t;
    }

    for (int i = 0; i < m; i++) {
        x_steps[i] = x;
        x += h;
    }

    return {
            t_steps, x_steps
    };
}

// =============================================================================================================
// =============================================================================================================
// KMB - Klasyczna metoda bezpośrednia
// =============================================================================================================
// =============================================================================================================

matrix_t get_KMB_solution(const int n, const int m) {
    matrix_t matrix = utils::matrix::get_matrix_with_conditions(n, m);

    for (int i = 1; i < n; i++) {
        for (int j = 1; j < m - 1; j++) {
            matrix[i][j] = matrix[i - 1][j] + KMB_LAMBDA *
                                              (matrix[i - 1][j - 1] - (2.0 * matrix[i - 1][j]) + matrix[i - 1][j + 1]);
        }
    }

    return matrix;
}

// =============================================================================================================
// =============================================================================================================
// Laasonen method using Thomas Algorithm
// =============================================================================================================
// =============================================================================================================


void calculate_eta(vector_t Lower, vector_t Diagonal, vector_t Upper, vector_t l, const int m) {
    for (int i = 1; i < m; i++) {
        l[i] = Lower[i] * (1 / Diagonal[i - 1]);
        Diagonal[i] -= l[i] * Upper[i-1];
    }
}

void calculate_r(vector_t b, vector_t l, const int m) {
    for (int i = 1; i < m; i++) {
        b[i] -= l[i] * b[i - 1];
    }
}

void calculate_x(vector_t Diagonal, vector_t Upper, vector_t b, vector_t x, const int m) {
    x[m - 1] = (1 / Diagonal[m - 1] * b[m - 1]);
    for (int i = m - 2; i >= 0; i--) {
        x[i] = (1 / Diagonal[i]) * (b[i] - Upper[i] * x[i + 1]);
    }
}


void Thomas(vector_t Lower, vector_t Diagonal, vector_t Upper, vector_t b, vector_t x, const int m) {
    vector_t vector_l = new double[m];
    calculate_eta(Lower, Diagonal, Upper, vector_l, m);
    calculate_r(b, vector_l, m);
    calculate_x(Diagonal, Upper, b, x, m);
    delete[] vector_l;
}

matrix_t get_Laasonen_Thomas(const int n, const int m) {
    matrix_t matrix_A = utils::matrix::get_matrix_with_conditions(n, m);

    const double lambda = 1.0 + 2.0 * LAASONEN_LAMBDA;

    vector_t Lower    = new double[m];
    vector_t Diagonal = new double[m];
    vector_t Upper    = new double[m];
    vector_t vector_b = new double[m];
    vector_t vector_x = new double[m];

    for (int k = 1; k < n; k++) {
        Lower[0]        = 0.0;
        Diagonal[0]     = 1.0;
        Upper[0]        = 0.0;
        vector_b[0]     = matrix_A[k - 1][0];

        for (int i = 1; i < m - 1; i++) {
            Lower[i]    = LAASONEN_LAMBDA;
            Diagonal[i] = -lambda;
            Upper[i]    = LAASONEN_LAMBDA;
            vector_b[i] = -matrix_A[k - 1][i];
        }

        Lower[m - 1]    = 0.0;
        Diagonal[m - 1] = 1.0;
        Upper[m - 1]    = 0.0;
        vector_b[m - 1] = matrix_A[k - 1][m - 1];


        Thomas(Lower, Diagonal, Upper, vector_b, vector_x, m);

        for (int i = 1; i < m; i++) {
            matrix_A[k][i] = vector_x[i];
        }
    }
    return matrix_A;
}

// =============================================================================================================
// =============================================================================================================
// Laasonen method using LU decomposition
// =============================================================================================================
// =============================================================================================================

int* get_vector_of_indices(const int m) {
    int* vec = new int[m];

    for (int i = 0; i < m; i++) {
        vec[i] = i;
    }

    return vec;
}

void swap_rows_if_necessary(matrix_t matrix, int k, const int m, int* vec_of_indices) {
    if (matrix[k][k] == 0) {
        int current_max_row = k + 1;
        double current_max = fabs(matrix[k+1][0]);

        for (int index = k + 1; index < m; index++) {
            double current_value = matrix[index][k];
            if (fabs(current_value) > current_max) {
                current_max = current_value;
                current_max_row = index;
            }
        }

        std::swap(matrix[current_max_row],matrix[k]);
        std::swap(vec_of_indices[current_max_row], vec_of_indices[k]);
    }
}

void calculate_next_iteration(matrix_t matrix, int k, const int m) {
    for (int i = k; i < m - 1; i++) {
        double reducer = (matrix[i + 1][k] / matrix[k][k]);
        for (int j = k; j < m; j++) {
            matrix[i + 1][j] -= matrix[k][j] * reducer;
        }
        matrix[i + 1][k] = reducer;
    }
}

void calculate_y(matrix_t matrix, vector_t y, vector_t b, const int m, int* vec_of_indices) {
    for (int i = 0; i < m; i++) {

        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += y[j] * matrix[i][j];
        }

        y[i] = b[vec_of_indices[i]] - sum;
    }
}

void LU_calculate_x(matrix_t matrix, vector_t y, const int m, vector_t x) {
    for (int i = m - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = m - 1; j > i; j--) {
            sum += x[j] * matrix[i][j];
        }

        x[i] = (y[i] - sum) / matrix[i][i];
    }
}

void LU(vector_t Lower, vector_t Diagonal, vector_t Upper, vector_t b, vector_t x, const int m) {
    vector_t y = new double[m];
    for (int i = 0; i < m; i++) {
        y[i] = 0.0;
    }

    int* vec_of_indices = get_vector_of_indices(m);

    matrix_t matrix = new vector_t[m];
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[m];
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            matrix[i][j] = 0.0;
        }
    }

    // upper
    for (int i = 0; i < m - 1; i++) {
        matrix[i][i + 1] = Upper[i];
    }

    // diagonal
    for (int i = 0; i < m; i++) {
        matrix[i][i] = Diagonal[i];
    }

    // lower
    for (int i = 1; i < m; i++) {
        matrix[i][i - 1] = Lower[i];
    }

    for (int i = 0; i < m - 1; i++) {
        swap_rows_if_necessary(matrix, i, m, vec_of_indices);
        calculate_next_iteration(matrix, i, m);
    }

    calculate_y(matrix, y, b, m, vec_of_indices);
    LU_calculate_x(matrix, y, m, x);

    for (int i = 0; i < m; i++) {
        delete[] matrix[i];
    }

    delete[] y;
    delete[] vec_of_indices;
    delete[] matrix;
}

matrix_t get_Laasonen_LU(const int n, const int m) {
    matrix_t matrix_A = utils::matrix::get_matrix_with_conditions(n, m);

    const double lambda = 1.0 + 2.0 * LAASONEN_LAMBDA;

    vector_t Lower    = new double[m];
    vector_t Diagonal = new double[m];
    vector_t Upper    = new double[m];
    vector_t vector_b = new double[m];
    vector_t vector_x = new double[m];

    for (int i = 0; i < m; i++) {
        vector_x[i] = 0.0;
    }

    for (int k = 1; k < n; k++) {
        Lower[0]        = 0.0;
        Diagonal[0]     = 1.0;
        Upper[0]        = 0.0;
        vector_b[0]     = matrix_A[k - 1][0];

        for (int i = 1; i < m - 1; i++) {
            Lower[i]    = LAASONEN_LAMBDA;
            Diagonal[i] = -lambda;
            Upper[i]    = LAASONEN_LAMBDA;
            vector_b[i] = -matrix_A[k - 1][i];
        }

        Lower[m - 1]    = 0.0;
        Diagonal[m - 1] = 1.0;
        Upper[m - 1]    = 0.0;
        vector_b[m - 1] = matrix_A[k - 1][m - 1];


        LU(Lower, Diagonal, Upper, vector_b, vector_x, m);

        for (int i = 1; i < m; i++) {
            matrix_A[k][i] = vector_x[i];
        }
    }

    return matrix_A;
}

// =============================================================================================================
// =============================================================================================================
//  Printing stuff for lab requirements. Later i make graphs from this data.
// =============================================================================================================
// =============================================================================================================

void KMB() {
    const double delta_t = get_delta_t(KMB_LAMBDA, h, D);
    const int n = ((t_max - t_min) / delta_t) + 2;
    const int m = ((x_end - x_start) / h);
    matrix_t analytic = get_analytic_solution(n, m, h, delta_t);
    matrix_t results = get_KMB_solution(n, m);

    Errors err = get_error(analytic, results, n, m);
    Steps steps = get_steps(delta_t, n, m);

    utils::matrix::write(analytic, n, m, "kmb_analytic.txt");
    utils::matrix::write(results, n, m, "kmb_solution.txt");
    utils::matrix::write(err.errors_matrix, n, m, "kmb_errors_matrix.txt");
    utils::vector::write(err.max_error, n, "kmb_max_error.txt");
    utils::vector::write(steps.t_steps, n, "kmb_t_steps.txt");
    utils::vector::write(steps.x_steps, m, "kmb_x_steps.txt");

    utils::cleanup(analytic, results, err, steps, n);
}

void Laasonen_Thomas() {
    const double delta_t = get_delta_t(LAASONEN_LAMBDA, h, D);
    const int n = ((t_max - t_min) / delta_t) + 2;
    const int m = ((x_end - x_start) / h);
    matrix_t analytic = get_analytic_solution(n, m, h, delta_t);
    matrix_t results = get_Laasonen_Thomas(n, m);

    Errors err = get_error(analytic, results, n, m);
    Steps steps = get_steps(delta_t, n, m);

    utils::matrix::write(analytic, n, m, "thomas_analytic.txt");
    utils::matrix::write(results, n, m, "thomas_solution.txt");
    utils::matrix::write(err.errors_matrix, n, m, "thomas_errors_matrix.txt");
    utils::vector::write(err.max_error, n, "thomas_max_error.txt");
    utils::vector::write(steps.t_steps, n, "thomas_t_steps.txt");
    utils::vector::write(steps.x_steps, m, "thomas_x_steps.txt");

    utils::cleanup(analytic, results, err, steps, n);
}

void Laasonen_LU() {
    const double delta_t = get_delta_t(LAASONEN_LAMBDA, h, D);
    const int n = ((t_max - t_min) / delta_t) + 2;
    const int m = ((x_end - x_start) / h);
    matrix_t analytic = get_analytic_solution(n, m, h, delta_t);
    matrix_t results = get_Laasonen_LU(n, m);

    Errors err = get_error(analytic, results, n, m);
    Steps steps = get_steps(delta_t, n, m);

    utils::matrix::write(analytic, n, m, "LU_analytic.txt");
    utils::matrix::write(results, n, m, "LU_solution.txt");
    utils::matrix::write(err.errors_matrix, n, m, "LU_errors_matrix.txt");
    utils::vector::write(err.max_error, n, "LU_max_error.txt");
    utils::vector::write(steps.t_steps, n, "LU_t_steps.txt");
    utils::vector::write(steps.x_steps, m, "LU_x_steps.txt");

    utils::cleanup(analytic, results, err, steps, n);
}

int main() {
    KMB();
    Laasonen_Thomas();
    Laasonen_LU();
    return 0;
}
