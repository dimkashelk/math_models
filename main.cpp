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
    if (is_equal(A, 0.0))
    {
    }
  }
  return res;
}
int main()
{
  std::cout << "Hello, World!" << std::endl;
  return 0;
}
