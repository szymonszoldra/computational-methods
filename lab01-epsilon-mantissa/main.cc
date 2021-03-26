#include <cmath>
#include <iostream>

int main() {
  float epsilon_float = 1.0f;
  float temp_epsilon_float = 1.0f;

  while ((1.f + temp_epsilon_float) != 1.f) {
    epsilon_float = temp_epsilon_float;
    temp_epsilon_float /= 2.0f;
  }

  double epsilon_double = 1.0;
  double temp_epsilon_double = 1.0;

  while ((1.0 + temp_epsilon_double) != 1.0) {
    epsilon_double = temp_epsilon_double;
    temp_epsilon_double /= 2.0;
  }

  int mantissa_float = -log2(epsilon_float);
  std::cout << "[FLOAT] Mantissa bits: " << mantissa_float << std::endl;
  std::cout << "[FLOAT] Epsilon: " << epsilon_float << std::endl;

  int mantissa_double = -log2(epsilon_double);
  std::cout << "[DOUBLE] Mantissa bits: " << mantissa_double << std::endl;
  std::cout << "[DOUBLE] Epsilon: " << epsilon_double << std::endl;

  return 0;
}