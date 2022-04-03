#include <math.h>
#include <stdlib.h>
#include <ctime>

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
const double v = 1.0 / 256.0;
const double m = 1.0;
const double l = 1.0;
const double k = 1.0;
const double dt = 1e-3;
const double y_init = 0.1 * l;
const double total_t = 2.0 * l / v;
const int total_steps = total_t / dt;
double y[1] = {0};
double p[1] = {0};
double t = 0;
double x = 2 * l;
double g = 0;
double E = 0;
double y_min = 0;
double y_max = 0;
double J = 0;

double y_array[total_steps] = {0};
double p_array[total_steps] = {0};
double x_array[total_steps] = {0};
double t_array[total_steps] = {0};
double J_array[total_steps] = {0};


double force(double *x, double v, double t) {
  double dist = 0;
  dist = sqrt((2.0 - v * t) * (2.0 - v * t) + x[0] * x[0]);
  return -k * (dist - l) * x[0] / dist;
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

  for (int i = 0; i < total_steps; i++) {
    t_array[i] = t;
    x_array[i] = x;
    y_array[i] = y[0];
    p_array[i] = p[0];
    E = energy(y[0], p[0]);
    if (E - potential(0.0) > 0) {
      y_min = dekker_root(tunnel_points, -2, 0, 1e-5);
      y_max = dekker_root(tunnel_points, 0, 2, 1e-5);
    }  // only two tunnelling points
    else {
      double y_M = 0;  // find the minimium of V to determine the root
      y_M = golden_section_opt(potential_inverse, 0, 2, 1e-5);
      y_min = dekker_root(tunnel_points, 0, y_M, 1e-5);
      y_max = dekker_root(tunnel_points, y_M, 2, 1e-5);
    }
    J = chebyshev_second_integral(J_func, y_min, y_max, 2e2);
    if (i % 5000 == 0) {
      cout << "steps" << i << ","
           << "total is " << total_steps << endl;
    }
    J_array[i] = J;
    Velocity_Verlet(y, p, force, dt, v, t);
    t += dt;
    x = x - v * dt;
  }
  stringstream fmt1;  // storage trajectories
  fmt1 << "2_trajectory_v=" << v << ".txt";
  ofstream OutFile1(fmt1.str());
  for (int i = 0; i < total_steps; i++) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << t_array[i] << " " << x_array[i] << " " << y_array[i] << " "
             << p_array[i] << " " << J_array[i] << endl;
  }
  OutFile1.close();
  fmt1.clear();
  return 0;
}