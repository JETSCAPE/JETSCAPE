#include <cmath>
#include <iostream>
#include <random>
#include <string>

// Generate Woods-Saxon numbers using the same method as in actual trento code.

int main(int /* argc */, char* argv[]) {
  double R;
  double a;
  int N;

  try {
    R = std::stod(std::string{argv[1]});
    a = std::stod(std::string{argv[2]});
    N = std::stoi(std::string{argv[3]});
  }
  catch (const std::exception&) {
    std::cerr << "usage: " << argv[0] << " R a n_samples\n";
    return 1;
  }

  std::mt19937_64 engine{std::random_device{}()};
  std::piecewise_linear_distribution<double>
    woods_saxon_dist{1000, 0, R+10.*a, [&R, &a](double r) {
      return r*r/(1. + std::exp((r-R)/a));
    }};

  for (int i = 0; i < N; ++i)
    std::cout << woods_saxon_dist(engine) << '\n';

  return 0;
}
