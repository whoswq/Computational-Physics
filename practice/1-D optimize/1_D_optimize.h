#ifndef _1_D_OPTIMAL_H_
#define _1_D_OPTIMAL_H_

const double PI = 3.1415926535897932384626433832;
double bisection_root_recursion(double (*func)(double), double a, double b,
                                double f_a, double f_b, double eps);
double bisection_root_loop(double (*func)(double), double a, double b,
                           double eps);
double newton_root(double (*func)(double), double (*dfunc)(double), double x,
                   double eps, int max_itr);
double secant_root(double (*func)(double), double x_0, double eps, int max_itr);
double dekker_root(double (*func)(double), double a, double b, double eps);
double brent_root(double (*func)(double), double a, double b, double eps);
double dekker_root_param(double (*func)(double *), double *var, int var_num,
                        double a, double b, double eps);
double golden_section_opt(double (*func)(double), double x_min, double x_max,
                          double eps);
double golden_section_opt_param(double (*func)(double *), double *var,
                                int var_num, double x_min, double x_max,
                                double eps);
double parabola_opt(double (*func)(double), double x1, double x2, double x3,
                    double eps);
double parabola_opt_param(double (*func)(double *), double *var, int var_num, double x1,
                    double x2, double x3, double eps);

#endif