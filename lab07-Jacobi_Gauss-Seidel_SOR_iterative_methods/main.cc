#include <iostream>
#include <array>
#include <cmath>

#define SIZE 4
#define TOLX 1e-7
#define TOLF 1e-7
#define MAX_ITER_COUNT 50

typedef std::array<std::array<double, SIZE>, SIZE> matrix_t;
typedef std::array<double, SIZE> vector_t;

matrix_t get_matrix_a() {
    return {{
        {100.0, -1.0, 2.0, -3.0},
        {1.0, 200.0, -4.0, 5.0},
        {-2.0, 4.0, 300.0, -6.0},
        {3.0, -5.0, 6.0, 400.0}
    }};
}

vector_t get_vector_b() {
    return {
        116.0,
        -226.0,
        912.0,
        -1174.0
    };
}

vector_t get_first_approximation() {
    return {
        2.0,
        2.0,
        2.0,
        2.0
    };
}

void print_vector(const vector_t &vector) {
    for (auto i: vector) {
        printf("%.9f ", i);
    }
}

bool check_and_print_en(const vector_t &xnplus1, const vector_t &xn) {
    int counter = 0;
    printf("\t Estymatory: [");
    for (int i = 0; i < SIZE; i++) {
        printf("%.9f ", fabs(xnplus1[i] - xn[i]));
        if (fabs(xnplus1[i] - xn[i]) < TOLX)
            counter++;
    }
    printf("]");
    return counter == SIZE;
}

bool check_and_print_residue(const matrix_t &A, const vector_t &b, const vector_t &xn) {
    int counter = 0;
    printf("\t Residuum: [");
    for (int i = 0; i < SIZE; i++) {
        double temp = 0;
        for (int j = 0; j < SIZE; j++) {
            temp += A[i][j] * xn[j];
        }
        printf("%.9f ", fabs(temp - b[i]));
        if (fabs(temp - b[i]) < TOLF)
            counter++;
    }
    printf("]\n");
    return counter == SIZE;
}

void jacobi_method() {
    std::cout << "Metoda Jacobiego" << std::endl;
    vector_t xn = get_first_approximation();
    vector_t b = get_vector_b();
    matrix_t A = get_matrix_a();

    vector_t xnplus1 = {0, 0, 0, 0};

    for (int iter_limit = 0; iter_limit < MAX_ITER_COUNT; iter_limit++) {
        for (int i = 0; i < SIZE; i++) {
            double temp = 0;

            for (int j = 0; j < SIZE; j++) {
                if (j == i)
                    continue;

                temp += A[i][j] * xn[j];
            }

            xnplus1[i] = (-temp + b[i]) / A[i][i];
        }

        printf("%2d    ", iter_limit);
        print_vector(xnplus1);
        bool en = check_and_print_en(xnplus1, xn);
        bool residue = check_and_print_residue(A, b, xn);
        if ( en && residue )
            break;
        std::swap(xnplus1, xn);
    }
}

void gauss_seidel_method() {
    std::cout << "Metoda Gaussa-Seidela" << std::endl;
    vector_t xn = get_first_approximation();
    vector_t b = get_vector_b();
    matrix_t A = get_matrix_a();
    vector_t right_side = {0, 0, 0, 0};

    for (int iter_limit = 0; iter_limit < MAX_ITER_COUNT; iter_limit++) {
        for (int i = 0; i < SIZE; i++) {
            double temp = 0;

            for (int j = i + 1 ; j < SIZE; j++) {
                temp += A[i][j] * xn[j];
            }

            right_side[i] = -temp + b[i];
        }

        vector_t xnplus1 = {0, 0, 0, 0};

        for (int i = 0; i < SIZE; i++) {
            double temp = 0;

            for (int j = 0; j <= i; j++) {
                temp += xnplus1[j] * A[i][j];
            }
            xnplus1[i] = (right_side[i] - temp) / A[i][i];
        }


        printf("%2d    ", iter_limit);
        print_vector(xnplus1);
        bool en = check_and_print_en(xnplus1, xn);
        bool residue = check_and_print_residue(A, b, xn);
        if ( en && residue )
            break;
        std::swap(xnplus1, xn);
    }
}

void sor_method() {
    std::cout << "Metoda SOR" << std::endl;
    vector_t xn = get_first_approximation();
    vector_t b = get_vector_b();
    matrix_t A = get_matrix_a();
    vector_t right_side = {0, 0, 0, 0};
    double OMEGA = 0.5;

    for (int iter_limit = 0; iter_limit < MAX_ITER_COUNT; iter_limit++) {
        for (int i = 0; i < SIZE; i++) {
            double temp = ((1.0 - 1.0 / OMEGA) * A[i][i]) * xn[i];

            for (int j = i + 1 ; j < SIZE; j++) {
                temp += A[i][j] * xn[j];
            }

            right_side[i] = -temp + b[i];
        }

        vector_t xnplus1 = {0, 0, 0, 0};

        for (int i = 0; i < SIZE; i++) {
            double temp = 0;

            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    temp += xnplus1[j] * (1.0 / OMEGA) * A[i][i];
                    continue;
                }
                temp += xnplus1[j] * A[i][j];
            }
            xnplus1[i] = (right_side[i] - temp) / ((1.0 / OMEGA) * A[i][i]);
        }


        printf("%2d    ", iter_limit);
        print_vector(xnplus1);
        bool en = check_and_print_en(xnplus1, xn);
        bool residue = check_and_print_residue(A, b, xn);
        if ( en && residue )
            break;
        std::swap(xnplus1, xn);
    }
}

int main() {
    jacobi_method();
    gauss_seidel_method();
    sor_method();
    return 0;
}
