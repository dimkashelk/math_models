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
double get_x2(double p1, double p2, double p4, double x1)
{
  // OK
  if (is_equal(p1, x1))
  {
    throw std::runtime_error("Check parameters");
  }
  return (p2 * p4 * x1 + x1 * x1 - x1) / (p1 - x1);
}
double get_x3(double p4, double x1)
{
  // OK
  return x1 / (1 + p4);
}
double get_p6(double p1, double p3, double p4, double p5, double x1, double x2, double x3)
{
  // OK
  return x2 - (-p1 * x2 - x1 * x2 + p5 * x3) / (p3 * p4);
}
void update_p4(double &p4)
{
  if (1.0 <= p4 && p4 < 2.0)
  {
    p4 += 0.1;
  }
  else if (2 <= p4 && p4 < 10.0)
  {
    p4 += 1;
  }
  else if (10 <= p4 && p4 < 100.0)
  {
    p4 += 10;
  }
  else if (100 <= p4 && p4 <= 1500.0)
  {
    p4 += 100;
  }
}
int main()
{
  double p1 = 8.4E-6, p2 = 6.6667E-4, p3 = 1.7778E-5, p5 = 2;
  std::vector<std::pair<double, double>> data_graphics;
  for (double p4 = 1.0; p4 <= 1500.0; update_p4(p4))
  {
    double a = -2 * p4 - 2;
    double b = 1 - p2 * p4 * p4 - 2 * p3 * p4 - 2 * p1 * p4 - p2 * p4 - 2 * p3 * p4 + p4 - 2 * p1 - p5;
    double c = -p3 * p4 * p4 - 2 * p1 * p4 - p3 * p4 - 2 * p1;
    double d = -p2 * p3 * p4 * p4 * p4 - p1 * p2 * p4 * p4 - p2 * p3 * p4 * p4 + p3 * p4 * p4 + p1 * p4 - p1 * p2 * p4 + p3 * p4 + p1 + p1 * p5;
    double first = -a;
    double second = a * p1 - b + c;
    double third = b * p1 + c * p2 * p4 - c - d;
    double fifth = d * p1;
    auto res = get_roots_cubic_equation(first, second, third, fifth);
    std::cout << "p4: " << p4 << ":\n";
    for (auto i: res)
    {
      double x2 = get_x2(p1, p2, p4, i);
      double x3 = get_x3(p4, i);
      std::cout << "\tx1: " << i << "\n";
      std::cout << "\t\tx2: " << x2 << "\n";
      std::cout << "\t\tx3: " << x3 << "\n";
      std::cout << "\t\tp6: " << get_p6(p1, p3, p4, p5, i, x2, x3) << "\n";
    }
  }
  return 0;
}
