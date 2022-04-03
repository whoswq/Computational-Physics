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
const double v = 1.0 / 1000.0;
const double m = 1.0;
const double l = 1.0;
const double k = 1.0;
const double dt = 1e-1;
const double y_init = -2.1;  // change this    
const double p_init = 0;
const double total_t = l / v;
const int total_steps = total_t / dt + 1;
const double delta = 1e-12;
const double amp = 2.0;  // change this
double y[1] = {0};
double p[1] = {0};
double t = 0;
double x = 0;
double g = amp * k * l / m;
double E = 0;
double y_min = 0;
double y_max = 0;
double J = 0;
double y_l[1] = {0};
double y_r[1] = {0};
double T = 0;

double y_array[total_steps] = {0};
double p_array[total_steps] = {0};
double g_array[total_steps] = {0};
double t_array[total_steps] = {0};
double J_array[total_steps] = {0};
double T_array[total_steps] = {0};
double x_array[total_steps] = {0};
double omega_array[total_steps] = {0};

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

double T_function(double y) {
  double up = 2 * m * (y - y_min) * (y_max - y);
  double down = E - potential(y);
  return sqrt(up / down);
}

double tunnel_points(double y) { return E - potential(y); }

void Velocity_Verlet(double *x, double *p, double (*f)(double *), double dt) {
  double p1 = 0;
  p1 = p[0] + 0.5 * dt * f(x);
  x[0] = p1 / m * dt + x[0];
  p[0] = p1 + 0.5 * dt * f(x);
}

int main() {
  y[0] = y_init;  // inintiate y pointer
  p[0] = p_init;
  clock_t start, end;
  start = clock();
  for (int i = 0; i < total_steps; i++) {
    t_array[i] = t;
    g_array[i] = g;
    y_array[i] = y[0];
    p_array[i] = p[0];
    x_array[i] = x;
    E = energy(y[0], p[0]);
    if (abs(x) <= l) {
      y_l[0] = -pow(abs(x) / l, 1.0 / 3.0) * sqrt(1.0 - pow(abs(x) / l, 2.0 / 3.0)) - delta;
      y_r[0] = -y_l[0];
      if (force(y_l) * force(y_r) > 0) {
        // only has one stable point
        double y_m = golden_section_opt(potential_inverse, -10.0, 10.0, 1e-7);
        y_min = dekker_root(tunnel_points, -10.0, y_m, 1e-10);
        y_max = dekker_root(tunnel_points, y_m, +10.0, 1e-10);
        T = chebyshev_anomaly_integral(T_function, y_min, y_max, 500);
        // the accuracy of chebyshev integral relies on the accuracy of root
        // especially the points increase
        J = chebyshev_second_integral(J_func, y_min, y_max, 500);
      } else if (force(y_l) * force(y_r) < 0) {
        // has two stable points
        double y_M =
            golden_section_opt(potential, y_l[0], y_r[0], 1e-7);
        double E_M = potential(y_M);
        if (E_M - E > delta) {
          // can not cross the barrier
          if (y[0] < y_M) {
            double y_m = golden_section_opt(potential_inverse, -10, y_M, 1e-7);
            y_min = dekker_root(tunnel_points, -10, y_m, 1e-10);
            y_max = dekker_root(tunnel_points, y_m, y_M, 1e-10);
            T = chebyshev_anomaly_integral(T_function, y_min, y_max, 500);
            J = chebyshev_second_integral(J_func, y_min, y_max, 500);
          } else if (y[0] > y_M) {
            double y_m = golden_section_opt(potential_inverse, y_M, 10, 1e-7);
            y_min = dekker_root(tunnel_points, y_M, y_m, 1e-10);
            y_max = dekker_root(tunnel_points, y_m, 10, 1e-10);
            T = chebyshev_anomaly_integral(T_function, y_min, y_max, 500);
            J = chebyshev_second_integral(J_func, y_min, y_max, 500);
          }
        } else if (E - E_M > delta) {
          y_min = dekker_root(tunnel_points, -10, y_M, 1e-10);
          y_max = dekker_root(tunnel_points, y_M, 10, 1e-10);
          T = chebyshev_anomaly_integral(T_function, y_min, y_max, 500);
          J = chebyshev_second_integral(J_func, y_min, y_max, 500);
        }
      }
    } else {
      // only has one stable point
      double y_m = golden_section_opt(potential_inverse, -10.0, 10.0, 1e-7);
      y_min = dekker_root(tunnel_points, -10.0, y_m, 1e-10);
      y_max = dekker_root(tunnel_points, y_m, +10.0, 1e-10);
      T = chebyshev_anomaly_integral(T_function, y_min, y_max, 500);
      J = chebyshev_second_integral(J_func, y_min, y_max, 500);
    }
    J_array[i] = J;
    omega_array[i] = 2.0 * PI / T;
    Velocity_Verlet(y, p, force, dt);
    t += dt;
    g = amp * k * l / m * cos(2.0 * PI * v * t);
    x = amp * l * sin(2.0 * PI * v * t);
    if (i % 50000 == 0) {
      end = clock();
      cout << "steps" << i << ","
           << "total is " << total_steps << ", time is "
           << (double)(end - start) / CLOCKS_PER_SEC << endl;
      cout << "J = " << J << endl;
      cout << "omega = " << 2.0 * PI / T << endl;
    }
  }
  stringstream fmt1;  // storage trajectories
  fmt1 << "4_berryphase_v=" << v << ".txt";
  ofstream OutFile1(fmt1.str());
  for (int i = 0; i < total_steps; i += 1) {
    OutFile1 << std::setiosflags(std::ios::scientific | std::ios::showpos)
             << t_array[i] << " " << g_array[i] << " " << x_array[i] << " "
             << y_array[i] << " " << p_array[i] << " " << J_array[i] << " "
             << omega_array[i] << endl;
  }
  OutFile1.close();
  fmt1.clear();
  return 0;
}