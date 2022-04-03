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
// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

double test_root_1(double x) { return sin(3.0 * x) + cos(7.0 * x); }
double dtest_root_1(double x) {
  return 3.0 * cos(2.0 * x) - 7.0 * sin(7.0 * x);
}
double test_root_2(double x) { return x * exp(-1.0 / (x * x)); }
double test_opt_1(double x) { return -x * x + 2.0 * x; }
double test_opt_2(double *var) {
  double x = var[0];
  double parm = var[1];
  return test_root_2(x) + parm * x;
}

double V_hat(double *var) {
  const double theta = var[0];
  const double lambda = var[1];
  const double mu = var[2];
  return ((1.0 - cos(theta)) - (1.0 - cos(theta) / mu) * (1.0 - cos(theta) / mu) /
                                 (2.0 * lambda * sin(theta) * sin(theta)));
}

int main() {
  const double lambda_min = 0;
  const double lambda_max = 1;
  const int lambda_num = 50;
  const double lambda_delta = (lambda_max - lambda_min) / ((double)lambda_num);
  const double mu = 0.5;
  double var[3] = {0};
  var[2] = mu;
  double lambda = 0;
  for (int i = 0; i < lambda_num; i++){
    lambda = lambda_min + (i + 0.5) * lambda_delta;
    var[1] = lambda;
    cout << golden_section_opt_param(V_hat, var, 3, 0, PI, 1e-6) << endl;
  } return 0;
}