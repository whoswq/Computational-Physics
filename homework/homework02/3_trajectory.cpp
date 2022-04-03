#include <math.h>
#include <stdlib.h>

#include <cmath>
#include <ctime>
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
const double v = 1.0 / 256.0;
const double m = 1.0;
const double l = 1.0;
const double k = 1.0;
const double dt = 2e-4;
const double y_init = -2 * l;
const double total_t = 0.5 * l / v;
const double delta = 1e-9;
const int total_steps = total_t / dt;
double y[1] = {0};
double p[1] = {0};
double t = 0;
double x = 0.2 * l;
double g = 2.0 * k * l / m;
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

double force(double *y) {
  double dist = 0;
  dist = sqrt(x * x + y[0] * y[0]);
  // to aviod divided by  zero
  if (y[0] > 0) {
    return -k * (dist - l) * (y[0] + delta) / (dist + delta) - m * g;
  } else {
    return -k * (dist - l) * (y[0] - delta) / (dist - delta) - m * g;
  }
}
double potential(double y) {
  double dis = sqrt(x * x + y * y) - l;
  return 0.5 * k * dis * dis + m * g * y;
}

double potential_inverse(double y) { return -potential(y); }

double J_func(double y) {
  double res =
      1.0 / PI * sqrt(2.0 * (E - potential(y)) / ((y - y_min) * (y_max - y)));
  return res;
}

double energy(double y, double p) { return p * p / (2.0 * m) + potential(y); }

double tunnel_points(double y) { return E - potential(y); }

void Velocity_Verlet(double *x, double *p, double (*f)(double *), double dt) {
  double p1 = 0;
  p1 = p[0] + 0.5 * dt * f(x);
  x[0] = p1 / m * dt + x[0];
  p[0] = p1 + 0.5 * dt * f(x);
}

int main() {
  y[0] = y_init;  // inintiate y pointer
  clock_t start, end;
  start = clock();
  for (int i = 0; i < total_steps; i++) {
    t_array[i] = t;
    g_array[i] = g;
    y_array[i] = y[0];
    p_array[i] = p[0];
    E = energy(y[0], p[0]);
    if (force(y_l) * force(y_r) > delta) {
      // only has one stable point
      double y_m = golden_section_opt(potential_inverse, -10.0, 10.0, 1e-4);
      y_min = dekker_root(tunnel_points, -10.0, y_m, 1e-6);
      y_max = dekker_root(tunnel_points, y_m, +10.0, 1e-6);
      J = chebyshev_second_integral(J_func, y_min, y_max, 1e3);
    } else if (force(y_l) * force(y_r) < -delta) {
      // has two stable points
      double y_M = golden_section_opt(potential, y_l[0], y_r[0], 1e-4);
      double E_M = potential(y_M);
      if (E_M - E > delta) {
        // can not cross the barrier
        if (y[0] < y_M) {
          double y_m = golden_section_opt(potential_inverse, -10, y_M, 1e-4);
          y_min = dekker_root(tunnel_points, -10, y_m, 1e-6);
          y_max = dekker_root(tunnel_points, y_m, y_M, 1e-6);
          J = chebyshev_second_integral(J_func, y_min, y_max, 1e3);
        } else if (y[0] > y_M) {
          double y_m = golden_section_opt(potential_inverse, y_M, 10, 1e-4);
          y_min = dekker_root(tunnel_points, y_M, y_m, 1e-6);
          y_max = dekker_root(tunnel_points, y_m, 10, 1e-6);
          J = chebyshev_second_integral(J_func, y_min, y_max, 1e3);
        }
      } else if (E - E_M > delta) {
        y_min = dekker_root(tunnel_points, -10, y_M, 1e-6);
        y_max = dekker_root(tunnel_points, y_M, 10, 1e-6);
        J = chebyshev_second_integral(J_func, y_min, y_max, 1e3);
      }
    }
    J_array[i] = J;
    Velocity_Verlet(y, p, force, dt);
    t += dt;
    g = 2.0 * k * l / m * cos(2 * PI * v * t);
    if (i % 5000 == 0) {
      end = clock();
      cout << "steps" << i << ","
           << "total is " << total_steps << ", time is "
           << (double)(end - start) / CLOCKS_PER_SEC << endl;
    }
  }
  stringstream fmt1;  // storage trajectories
  fmt1 << "3_trajectory_v=" << v << ".txt";
  ofstream OutFile1(fmt1.str());
  for (int i = 0; i < total_steps; i++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << t_array[i] << " " << g_array[i] << " " << y_array[i] << " "
             << p_array[i] << " " << J_array[i] << endl;
  }
  OutFile1.close();
  fmt1.clear();
  return 0;
}