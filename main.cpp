#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
bool is_equal(double x, double y)
{
  return std::fabs(x - y) < std::numeric_limits< double >::epsilon();
}
template < typename T >
int sign(T val)
{
  return (T(0) < val) - (val < T(0));
}
std::vector< double > get_roots_cubic_equation(double a, double b, double c, double d)
{
  if (is_equal(a, 0.0))
  {
    throw std::runtime_error("Check parameters");
  }
  b /= a;
  c /= a;
  d /= a;
  double Q = (std::pow(b, 2) - 3 * c) / 9;
  double R = (2 * std::pow(b, 3) - 9 * b * c + 27 * d) / 54;
  std::vector< double > res;
  if (std::pow(R, 2) < std::pow(Q, 3))
  {
    double t = std::acos(R / std::sqrt(std::pow(Q, 3))) / 3;
    res.push_back(-2 * std::sqrt(Q) * std::cos(t) - b / 3);
    res.push_back(-2 * std::sqrt(Q) * std::cos(t + (2 * M_PI / 3)) - b / 3);
    res.push_back(-2 * std::sqrt(Q) * std::cos(t - (2 * M_PI / 3)) - b / 3);
  }
  else
  {
    double A = -sign(R) * std::pow(std::fabs(R) + std::sqrt(std::pow(R, 2) - std::pow(Q, 3)), 1.0 / 3.0);
    double B = 0.0;
    if (!is_equal(A, 0.0))
    {
      B = Q / A;
    }
    res.push_back((A + B) - b / 3);
    if (is_equal(A, B))
    {
      res.push_back(-A - b / 3);
      res.push_back(-A - b / 3);
    }
  }
  return res;
}
double get_a(double p1, double p2, double p3, double p4)
{
  return std::pow(p1, 2) * p2 * p4 - p1 * p3 * p4 + p1 * p2 * p3 * std::pow(p4, 2) - p2 * p3 * std::pow(p4, 2);
}
double get_b(double p1, double p3, double p4)
{
  return 2 * p1 - p3 * p4 + 1;
}
double get_c(double p1, double p2, double p3, double p4)
{
  return 2 * std::pow(p1, 2) - 3 * p1 + 2 * p1 * p3 * p4 + 2 * p1 * p2 * p4 + p2 * p3 * std::pow(p4, 2) - 1;
}
double get_x2(double p1, double p2, double p4, double x1) {
  if (is_equal(p1, x1)) {
    throw std::runtime_error("Check parameters");
  }
  return (p2 * p4 * x1 + x1 * x1 - x1) / (p1 - x1);
}
int main()
{
  double p1 = 8.4E-6, p2 = 6.6667E-4, p3 = 1.7778E-5, p5 = 2;
  for (double p4 = 1.0; p4 <= 2.05; p4 += 0.1)
  {
    double parameter = (-1 - p4);
    double a = -2 * parameter;
    double b = get_b(p1, p3, p4) * parameter + p5;
    double c = get_c(p1, p2, p3, p4) * parameter - 2 * p1 * p5;
    double d = get_a(p1, p2, p3, p4) * parameter + p1 * p1 * p5;
    auto res = get_roots_cubic_equation(a, b, c, d);
    std::cout << p4 << ":\n";
    for (auto i: res) {
      std::cout << "\t" << i << "\n";
    }
  }
  return 0;
}
