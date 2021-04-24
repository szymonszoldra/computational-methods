#include <iostream>
#include <array>
#include <cmath>

#define SIZE 4

typedef std::array<std::array<double, SIZE>, SIZE> matrix_t;
typedef std::array<double, SIZE> vector_t;
typedef std::array<int, SIZE> index_vector_t;

index_vector_t vec_of_indexes = {0, 1, 2, 3};

matrix_t get_matrix_a() {
    matrix_t m = {{
            { 1.0, -20.0, 30.0, -4.0 },
            { 2.0, -40.0, -6.0, 50.0 },
            { 9.0, -180.0, 11.0, -12.0 },
            { -16.0, 15.0, -140.0, 13.0 }
    }};

    return m;
};

vector_t get_vector_b() {
    return { 35.0, 104.0, -366.0, -354.0 };
}

void print_matrix(const matrix_t &matrix, const char* name) {
    printf("\n %s \n", name);
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < 4; j++) {
            printf("%f\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void swap_rows_if_necessary(matrix_t &matrix, int k ) {
    if (matrix[k][k] == 0) {
        int current_max_row = k + 1;
        double current_max = fabs(matrix[k+1][0]);

        for (int index = k + 1; index < SIZE; index++) {
            double current_value = matrix[index][k];
            if (fabs(current_value) > current_max) {
                current_max = current_value;
                current_max_row = index;
            }
        }

        std::swap(matrix[current_max_row],matrix[k]);
        std::swap(vec_of_indexes[current_max_row], vec_of_indexes[k]);
    }
}

void calculate_next_iteration(matrix_t &matrix, int k) {
    for (int i = k; i < SIZE - 1; i++) {
        double reducer = (matrix[i + 1][k] / matrix[k][k]);
        for (int j = k; j < SIZE; j++) {
            matrix[i + 1][j] -= matrix[k][j] * reducer;
        }
        matrix[i + 1][k] = reducer;
    }
}

void LU_decomposition(matrix_t &matrix) {

    for (int k = 0; k < SIZE - 1; k++) {
        swap_rows_if_necessary(matrix, k);
        calculate_next_iteration(matrix, k);
    }

    print_matrix(matrix, "MATRIX");
}

vector_t calculate_y(matrix_t &matrix, vector_t &b) {
    vector_t y = {0, 0, 0, 0};

    for (int i = 0; i < SIZE; i++) {

        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += y[j] * matrix[i][j];
        }

        y[i] = b[vec_of_indexes[i]] - sum;
    }
    return y;
}

vector_t calculate_x(matrix_t &matrix, vector_t &y) {
    vector_t x = {0, 0, 0, 0};

    for (int i = SIZE - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = SIZE - 1; j > i; j--) {
            sum += x[j] * matrix[i][j];
        }

        x[i] = (y[i] - sum) / matrix[i][i];
    }

    return x;
}

void check(const matrix_t &A, const vector_t &b, const vector_t &x) {
    std::cout << "\nJesli ponizej są cztery zera to znaczy że wektor x został wyliczony poprawnie:" << std::endl;
    for (int i = 0; i < SIZE; i++) {
        int sum = 0;
        for (int j = 0; j < SIZE; j++) {
            sum += A[i][j] * x[j];
        }
        sum -= b[i];
        std::cout << sum << std::endl;
    }
}

int main() {

    matrix_t matrix_a = get_matrix_a();
    vector_t vector_b = get_vector_b();

    // This version of LU decomposition deos not create two separate arrays to hold the values of the Lower and Upper matrices.
    // Instead, I use the previously created matrix A to combine two matrices into one (there are Upper values on the main diagonal).
    // Later, when counting the vector Y, I take this fact into account, I calculate it as if there were 1 on the diagonal.
    
    LU_decomposition(matrix_a);

    vector_t vector_y = calculate_y(matrix_a, vector_b);
    std::cout << "Wektor Y:" << std::endl;
    for (int i = 0; i < SIZE; i++) {
        std::cout << "y" << i+1 << " = " << vector_y[i] << std::endl;
    }

    std::cout << std::endl;

    vector_t vector_x = calculate_x(matrix_a, vector_y);
    std::cout << "Wektor X:" << std::endl;
    for (int i = 0; i < SIZE; i++) {
        std::cout << "x" << i+1 << " = " << vector_x[i] << std::endl;
    }

    check(get_matrix_a(), get_vector_b(), vector_x);
    return 0;
}
