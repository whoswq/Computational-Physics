#include <math.h>
#include <stdlib.h>
#include <windows.h>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "1_D_optimize.h"
#include "integral.h"
// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

// need to change v
const double v = 1.0 / 64.0;
const double m = 1.0;
const double l = 1.0;
const double k = 1.0;
const double dt = 2e-4;
const double y_init = -1.0;
const double p_init = 2.5;
const double total_t = 0.5 * l / v;
const int total_steps = total_t / dt;
double y[1] = {0};
double p[1] = {0};
double t = 0;
double x = 0.2 * l;
double g = -2.0 * k * l / m;
double E = 0;
double y_min = 0;
double y_max = 0;
double J = 0;
double y_l[1] = {-x * sqrt(pow(l / x, 2.0 / 3.0) - 1.0)};
double y_r[1] = {-y_l[0]};

double y_array[total_steps] = {0};
double p_array[total_steps] = {0};
double g_array[total_steps] = {0};
double t_array[total_steps] = {0};
double J_array[total_steps] = {0};

double dekker_root(double (*func)(double), double a, double b, double eps) {
  double f_a = func(a);
  double f_b = func(b);
  double var = 0;
  double f_var = 0;
  double m = 0;
  double s = 0;
  double b_old = b - (b - a) / 4;  // generate a b' to apply secant method
  double fb_old = func(b_old);
  if (f_a * f_b > 0) {
    throw "error in dekker_root, no root";
    return 0.0;
  }
  while (abs(b - a) > eps and f_b != 0) {
    if (abs(f_a) < abs(f_b)) {
      // keep b is a better aproximation
      var = b;
      b = a;
      a = var;
      f_var = f_b;
      f_b = f_a;
      f_a = f_var;
    }
    m = (a + b) / 2.0;                                // mid=point
    s = (b * fb_old - b_old * f_b) / (fb_old - f_b);  // secant point
    b_old = b;
    // determine bisection or secant method
    if ((s - m) * (s - b) < 0) {  // secant point is between mid-point and b
      b = s;                      // choose secant point
    } else {
      b = m;  // choose mid-point
    }
    // update a, b
    fb_old = f_b;
    f_b = func(b);
    if (f_b * fb_old < 0) {  // root is between b and b_old
      a = b_old;
      f_a = fb_old;
    }
  }
  return b;
}

double parabola_opt(double (*func)(double), double x1, double x2, double x3,
                    double eps) {
  double f1 = func(x1);
  double f2 = func(x2);
  double f3 = func(x3);
  double x0 = 0;
  double f0 = 0;
  while ((abs(x1 - x2) > eps) and (abs(x1 - x3) > eps) and
         (abs(x2 - x3) > eps)) {
    x0 = ((f1 * (x2 * x2 - x3 * x3) + f2 * (x3 * x3 - x1 * x1) +
           f3 * (x1 * x1 - x2 * x2)) /
          (f1 * (x2 - x3) + f2 * (x3 - x1) + f3 * (x1 - x2))) /
         2.0;
    f0 = func(x0);
    if ((x0 - x1) * (x0 - x2) < 0) {  // x0 is between x1 & x2
      if (f0 > f2) {                  // drop x3
        x3 = x2;
        f3 = f2;
        x2 = x0;
        f2 = f0;
      } else {  // drop x1
        x1 = x0;
        f1 = f0;
      }
    } else {          // x0 is between x2 & x3
      if (f0 > f2) {  // drop x1
        x1 = x2;
        f1 = f2;
        x2 = x0;
        f2 = f0;
      } else {  // drop x3
        x3 = x0;
        f3 = f0;
      }
    }
  }
  if ((f1 > f2) and (f1 > f3)) {
    x0 = x1;
  } else if (f2 > f3) {
    x0 = x2;
  } else {
    x0 = x3;
  }
  return x0;
}

double golden_section_opt(double (*func)(double), double x_min, double x_max,
                          double eps) {
  double tau = (sqrt(5.0) - 1.0) / 2.0;
  double a = x_min;
  double b = x_max;
  double alpha = a + (b - a) * tau;
  double beta = a + b - alpha;
  double fa = func(alpha);
  double fb = func(beta);
  while (abs(b - a) > eps) {
    if (fb < fa) {
      a = beta;
      beta = alpha;
      fb = fa;                    // use previous result
      alpha = a + (b - a) * tau;  // new point
      fa = func(alpha);
    } else {
      b = alpha;
      alpha = beta;
      fa = fb;
      beta = b - (b - a) * tau;
      fb = func(beta);
    }
  }
  double res = 0;
  res = parabola_opt(func, a, alpha, b, eps / 5.0);
  return res;
}

double chebyshev_second_integral(double (*func)(double), double a, double b,
                                 int n) {
  double integral = 0;
  double x_i = 0;
  double cs = 0;
  for (int i = 1; i < n + 1; i++) {
    cs = cos(PI * (2 * i - 1) / (2 * n));
    x_i = (a + b) / 2 + (b - a) / 2 * cs;
    double f_i = func(x_i);
    integral += func(x_i) * (1 - cs * cs);
  }
  return PI / n * integral * (b - a) * (b - a) * 0.25;
}

double force(double *x, double v, double t) {
  double dist = 0;
  dist = sqrt((2.0 - v * t) * (2.0 - v * t) + x[0] * x[0]);
  return -k * (dist - l) * x[0] / dist - m * g;
}

double potential(double y) {
  double dis = sqrt(x * x + y * y) - l;
  return 0.5 * k * dis * dis + m * g * y;
}

double potential_inverse(double y) { return -potential(y); }

double J_func(double y) {
  double p = potential(y);
  double d = (y - y_min) * (y_max - y);
  return sqrt(2.0 * (E - potential(y)) / (y - y_min) * (y_max - y));
}

double energy(double y, double p) { return p * p / (2.0 * m) + potential(y); }

double tunnel_points(double y) { return E - potential(y); }

void Velocity_Verlet(double *x, double *p,
                     double (*f)(double *, double, double), double dt, double v,
                     double t) {
  double p1 = 0;
  p1 = p[0] + 0.5 * dt * f(x, v, t);
  x[0] = p1 / m * dt + x[0];
  p[0] = p1 + 0.5 * dt * f(x, v, t);
}

int main() {
  y[0] = y_init;  // inintiate y pointer
  p[0] = p_init;
  for (int i = 0; i < total_steps; i++) {
    t_array[i] = t;
    g_array[i] = g;
    y_array[i] = y[0];
    p_array[i] = p[0];
    E = energy(y[0], p[0]);
    Velocity_Verlet(y, p, force, dt, v, t);
    t += dt;
  }
  stringstream fmt1;  // storage trajectories
  fmt1 << "test" << ".txt";
  ofstream OutFile1(fmt1.str());
  for (int i = 0; i < total_steps; i++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << t_array[i] << " " << y_array[i] << " "
             << p_array[i] << endl;
  }
  OutFile1.close();
  fmt1.clear();
  return 0;
}