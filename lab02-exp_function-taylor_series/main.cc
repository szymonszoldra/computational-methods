#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

double get_relative_error(double approx_value, double value) {
  return fabs(((approx_value - value) / value));
}

double f(double x) { return ((1. - exp(-x)) / x); }

double taylor_series(double x) {
  double current_value = 1.0;
  double term = 1.0;

  for (int n = 2; n < 19; n++) {
    term *= -x / n;
    current_value += term;
  }

  return current_value;
}

int main() {
  std::fstream file("data.txt");
  std::ofstream output("error_plot.txt");
  std::vector<double> log10x;
  std::vector<double> x_values;
  std::vector<double> correct_answers;
  double input;

  int counter = 0;

  while (file >> input) {
    int modulo = counter % 3;
    switch (modulo) {
      case 0:
        log10x.push_back(input);
        break;
      case 1:
        x_values.push_back(input);
        break;
      case 2:
        correct_answers.push_back(input);
        break;
      default:
        break;
    }
    counter++;
  }

  for (int i = 0; i < counter / 3; i++) {
    if (log10x[i] < 0) {
      double answer = taylor_series(x_values[i]);
      double error = get_relative_error(answer, correct_answers[i]);
      double log_error = log10(error);
      output << log10x[i] << " " << log_error << std::endl;
    } else {
      double answer = f(x_values[i]);
      double error = get_relative_error(answer, correct_answers[i]);
      double log_error = log10(error);
      output << log10x[i] << " " << log_error << std::endl;
    }
  }
  return 0;
}