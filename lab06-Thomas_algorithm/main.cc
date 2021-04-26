#include <iostream>
#include <array>

#define SIZE 6

typedef std::array<std::array<double, 3>, SIZE> matrix_t;
typedef std::array<double, SIZE> vector_t;

vector_t vector_l = {0, 0, 0, 0, 0, 0};

matrix_t get_matrix_a() {
    return {{
                    {0, 10.0, 1.0 / 2.0},
                    {1.0 / 3.0, 20.0, 1.0 / 4.0},
                    {1.0 / 5.0, 30.0, 1.0 / 6.0},
                    {1.0 / 7.0, 30.0, 1.0 / 8.0},
                    {1.0 / 9.0, 20.0, 1.0 / 10.0},
                    {1.0 / 11.0, 10.0, 0},
            }};
}

vector_t get_vector_b() {
    return { 31.0,
             165.0 / 4.0,
             917.0 / 30.0,
             851.0 / 28.0,
             3637.0 / 90.0,
             332.0 / 11.0
            };
}

void calculate_eta(matrix_t &matrix) {
    for (int i = 1; i < SIZE; i++) {
        vector_l[i] = matrix[i][0] * (1 / matrix[i - 1][1]);
        matrix[i][1] -= vector_l[i] * matrix[i-1][2];
    }
}

void calculate_r(vector_t &vector_b) {
    for (int i = 1; i < SIZE; i++) {
        vector_b[i] -= vector_l[i] * vector_b[i - 1];
    }
}

void calculate_x(const matrix_t &matrix, const vector_t &vector_b, vector_t &vector_x) {
    for (int i = SIZE - 1; i >=0; i--) {
        vector_x[i] = (1 / matrix[i][1]) * (
                                            i == (SIZE - 1)
                                            ?  vector_b[i]
                                            : (vector_b[i] - matrix[i][2] * vector_x[i + 1])
                                            );
    }
}

void print_result(const vector_t &vector_x) {
    for (const double value: vector_x) {
        std::cout << value << std::endl;
    }
}

void thomas_algorithm(matrix_t &matrix, vector_t &vector_b, vector_t &vector_x) {
    calculate_eta(matrix);
    calculate_r(vector_b);
    calculate_x(matrix, vector_b, vector_x);
    print_result(vector_x);
}

int main() {
    matrix_t matrix = get_matrix_a();
    vector_t vector_b = get_vector_b();
    vector_t vector_x = {0, 0, 0, 0, 0, 0};

    thomas_algorithm(matrix, vector_b, vector_x);

    return 0;
}
