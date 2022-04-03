#include "1_D_optimize.h"

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
// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

double bisection_root_recursion(double (*func)(double), double a, double b,
                                double f_a, double f_b, double eps) {
  if (f_a * f_b > 0) {
    throw "error in bisection_root_recursion, no root";
    return 0.0;
  }
  double mid = (b + a) / 2;
  if (abs(b - a) < eps) {
    return mid;
  }
  double f_mid = func(mid);
  double root = 0.0;
  if (f_mid * f_a < 0) {
    root = bisection_root_recursion(func, a, mid, f_a, f_mid, eps);
  } else {
    root = bisection_root_recursion(func, mid, b, f_mid, f_b, eps);
  }
  return root;
}

double bisection_root_loop(double (*func)(double), double a, double b,
                           double eps) {
  double f_a = func(a);
  double f_b = func(b);
  if (f_a * f_b > 0) {
    throw "error in bisection_root_recursion, no root";
    return 0.0;
  }
  double mid = (a + b) / 2;
  double f_mid = func(mid);
  while (abs(a - b) > eps) {
    if (f_mid * f_a < 0) {
      f_b = f_mid;
      b = mid;
    } else {
      f_a = f_mid;
      a = mid;
    }
    mid = (a + b) / 2;
    f_mid = func(mid);
  }
  return mid;
}

double newton_root(double (*func)(double), double (*dfunc)(double), double x,
                   double eps, int max_itr = 1e5) {
  double f = func(x);
  int cnt = 0;
  while (abs(f) > eps and cnt < max_itr) {
    x = x - f / dfunc(x);
    f = func(x);
    cnt += 1;
  }
  if (cnt == max_itr - 1) {
    throw "error in newton_root, does not converge";
  }
  return x;
}

double secant_root(double (*func)(double), double x_0, double eps,
                   int max_itr = 1e5) {
  double x_pre = x_0;
  double x_nxt = x_0 / 2.0;
  double f_pre = func(x_pre);
  double f_nxt = func(x_nxt);
  double var = 0;
  int cnt = 0;

  while (abs(f_nxt) > eps) {
    var = x_nxt;
    x_nxt = (x_pre * f_nxt - x_nxt * f_pre) / (f_nxt - f_pre);
    x_pre = var;
    f_pre = func(x_pre);
    f_nxt = func(x_nxt);
  }
  if (cnt == max_itr - 1) {
    throw "error in secant_root, does not converge";
  }
  return x_nxt;
}

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

/*find root using dekker method, func can have parameters
func: return double, input double pointer
var: the input of func
var_num: the number of var
a: left end
b: right end
eps:
*/
double dekker_root_param(double (*func)(double *), double *var, int var_num,
                         double a, double b, double eps) {
  double var_copy[var_num] = {0.0};
  for (int i = 0; i < var_num; i++) {
    var_copy[i] = var[i];
  }
  var_copy[0] = a;
  double f_a = func(var_copy);
  var_copy[0] = b;
  double f_b = func(var_copy);
  double var_, f_var, m, s = 0;
  double b_old = b - (b - a) / 4;  // generate a b' to apply secant method
  var_copy[0] = b_old;
  double fb_old = func(var_copy);
  if (f_a * f_b > 0) {
    throw "error in dekker_root, no root";
    return 0.0;
  }
  while (abs(b - a) > eps and f_b != 0) {
    if (abs(f_a) < abs(f_b)) {
      // keep b is a better aproximation
      var_ = b;
      b = a;
      a = var_;
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
    var_copy[0] = b;
    f_b = func(var_copy);
    if (f_b * fb_old < 0) {  // root is between b and b_old
      a = b_old;
      f_a = fb_old;
    }
  }
  return b;
}

double brent_root(double (*func)(double), double a, double b, double eps) {
  double fa = func(a);
  double fb = func(b);
  double var, f_var, c, d, fc, s, fs = 0;
  bool mflag = true;
  if (fa * fb > 0) {
    throw "error in brent_root, no root";
    return 0.0;
  }
  if (abs(fa) < abs(fb)) {
    var = b;
    b = a;
    a = var;
    f_var = fb;
    fb = fa;
    fa = f_var;
  }
  c = a;
  d = a;
  fc = fa;
  mflag = true;
  while (abs(a - b) > eps and fb != 0) {
    if (fa != fc and fb != fc) {  // inverse quadratic interpolation
      s = a * fb * fc / ((fa - fb) * (fa - fc)) +
          b * fa * fc / ((fb - fa) * (fb - fc)) +
          c * fa * fb / ((fc - fa) * (fc - fb));
    } else {
      s = (a * fb - b * fa) / (fb - fa);
    }
    if ((s - (3.0 * a + b) / 4.0) * (s - b) > 0 or
        (mflag and (abs(s - b) >= abs(b - c) / 2.0)) or
        (!mflag and (abs(s - b) >= abs(c - d) / 2.0)) or
        (mflag and (abs(b - c) < eps)) or (!mflag and (abs(c - d) < eps))) {
      s = (a + b) / 2.0;
      mflag = true;
    } else {
      mflag = false;
    }
    fs = func(s);
    d = c;
    c = b;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    if (abs(fb) > abs(fa)) {
      var = b;
      f_var = fb;
      b = a;
      fb = fa;
      a = var;
      fa = f_var;
    }
  }
  return b;
}

// it can not guarantee convergency
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

double parabola_opt_param(double (*func)(double *), double *var, int var_num,
                          double x1, double x2, double x3, double eps) {
  double var_copy[var_num] = {0};
  for (int i = 0; i < var_num; i++) {
    var_copy[i] = var[i];
  }
  var_copy[0] = x1;
  double f1 = func(var_copy);
  var_copy[0] = x2;
  double f2 = func(var_copy);
  var_copy[0] = x3;
  double f3 = func(var_copy);
  double x0 = 0;
  double f0 = 0;
  while ((abs(x1 - x2) > eps) and (abs(x1 - x3) > eps) and
         (abs(x2 - x3) > eps)) {
    x0 = ((f1 * (x2 * x2 - x3 * x3) + f2 * (x3 * x3 - x1 * x1) +
           f3 * (x1 * x1 - x2 * x2)) /
          (f1 * (x2 - x3) + f2 * (x3 - x1) + f3 * (x1 - x2))) /
         2.0;
    var_copy[0] = x0;
    f0 = func(var_copy);
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
  res = parabola_opt(func, a, alpha, b, eps / 2.0);
  return res;
}

/*
func: return double, input double pionter
var: parameters of func
*/
double golden_section_opt_param(double (*func)(double *), double *var,
                                int var_num, double x_min, double x_max,
                                double eps) {
  double var_copy[var_num] = {0};
  for (int i = 0; i < var_num; i++) {
    var_copy[i] = var[i];
  }
  double tau = (sqrt(5.0) - 1.0) / 2.0;
  double a = x_min;
  double b = x_max;
  double alpha = a + (b - a) * tau;
  double beta = a + b - alpha;
  var_copy[0] = alpha;
  double fa = func(var_copy);
  var_copy[0] = beta;
  double fb = func(var_copy);
  while (abs(b - a) > eps) {
    if (fb < fa) {
      a = beta;
      beta = alpha;
      fb = fa;                    // use previous result
      alpha = a + (b - a) * tau;  // new point
      var_copy[0] = alpha;
      fa = func(var_copy);
    } else {
      b = alpha;
      alpha = beta;
      fa = fb;
      beta = b - (b - a) * tau;
      var_copy[0] = beta;
      fb = func(var_copy);
    }
  }
  double res = 0;
  res = parabola_opt_param(func, var, var_num, a, alpha, b, eps / 5.0);
  return res;
}

/*
double test_opt_1(double x) { return -x * x + 2.0 * x; }
double test_opt_2(double *var) {
  double x = var[0];
  double parm = var[1];
  return test_opt_1(x) + parm * x;
}

double V_hat(double *var) {
  const double theta = var[0];
  const double lambda = var[1];
  const double mu = var[2];
  return ((1.0 - cos(theta)) - (1.0 - cos(theta) / mu) *
                                   (1.0 - cos(theta) / mu) /
                                   (2.0 * lambda * sin(theta) * sin(theta)));
}

int main() {
  const double lambda = 1;
  const double mu = 2;
  double var[3] = {0};
  var[2] = mu;
  var[1] = lambda;
  cout << golden_section_opt_param(V_hat, var, 3, 0, PI, 1e-6) << endl;
  //cout << golden_section_opt_param(test_opt_2, var, 2, -3, 3, 1e-6) << endl;

  return 0;
}
*/
