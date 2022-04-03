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

const double m = 1.0;
const double l = 1.0;
const double k = 1.0;
const double dt = 1e-4;
const double y_init =
    -2.1;  // -2.1 * l for problem 4; 0.52 * l for problem 5
const double p_init = 0;
const double y_fina = -3.118229e+00;  // +8.542260e-01 -9.279862e-02 for problem
                                      // 5 -2.265301e+00 -5.198193e-01
                                      // -3.118229e+00 -8.912021e-01
const double p_fina =
    -8.912021e-01;     // -3.833860e+00 -3.387657e-01 for problem 4
const double amp = 2;  // 2 for problem 4; 0.3 for problem 5
const int total_steps = 10.0 / dt;
const double delta = 1e-12;
const double eps = 1e-3;  // 1e-4 for problem 5,
double y[1] = {0};
double p[1] = {0};
double t = 0;
double x = 0;
double g = amp * k * l / m;

double y_array[total_steps] = {0};
double p_array[total_steps] = {0};

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
    Velocity_Verlet(y, p, force, dt);
    t += dt;
    if (i % 50000 == 0) {
      end = clock();
      cout << "steps" << i << ","
           << "total is " << total_steps << ", time is "
           << (double)(end - start) / CLOCKS_PER_SEC << endl;
    }
    if (abs(y[0] - y_fina) < eps and abs(p[0] - p_fina) < eps) {
      cout << "total time is :" << t << endl;
    }
  }
  return 0;
}