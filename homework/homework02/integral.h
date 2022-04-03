#ifndef _INTEGRAL_H_
#define _INTEGRAL_H_

const double PI = 3.1415926535897932384626433832;

double midpoint_intergral(double (*func)(double), double a, double b, int n);
double trapezoid_integral(double (*func)(double), double a, double b, int n);
double self_adaptive_integral(double (*func)(double), double a, double b, int n, double err);
double refine(double (*func)(double), double fa, double fb, double a, double b, double err);
double simpson_integral(double (*func)(double), double a, double b, int n);
double richardson_integral(double (*func)(double), double a, double b, int n, int k);
double chebyshev_anomaly_integral(double (*func)(double), double a, double b, int n);
double chebyshev_second_integral(double (*func)(double), double a, double b, int n);
#endif