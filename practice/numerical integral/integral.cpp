#include "integral.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <random>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <windows.h>
// using namespace std;
using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

double midpoint_integral(double (*func)(double), double a, double b, int n)
{
    double dx = (b - a) / n;
    double integral = 0;
    double x = 0;
    for (int i = 0; i < n; i++)
    {
        x = ((double)i + 0.5) * dx + a;
        integral += dx * func(x);
    }
    return integral;
}

double trapezoid_integral(double (*func)(double), double a, double b, int n)
{
    double dx = (b - a) / n;
    double integral = (func(a) + func(b)) / 2 * dx;
    double x = 0;
    for (int i = 1; i < n; i++)
    {
        x = a + dx * i;
        integral += dx * func(x);
    }
    return integral;
}

double refine(double (*func)(double), double fa, double fb, double a, double b, double err = 1e-4)
{
    double c = (a + b) / 2;
    double fc = func(c);
    double integral = 0;
    if (abs(2.0 * fc - fa - fb) < err)
    {
        integral = (b - a) / 4 * (fa + 2.0 * fc + fb);
    }
    else
    {
        integral = refine(func, fa, fc, a, c, err) + refine(func, fc, fb, c, b, err);
    }
    return integral;
}

double self_adaptive_integral(double (*func)(double), double a, double b, int n, double err = 1e-4)
{
    double integral = 0;
    double dx = (b - a) / n;
    double f0 = func(a);
    double x0 = a;
    double x1 = 0;
    double f1 = 0;
    for (int i = 1; i < n + 1; i++)
    {
        x1 = a + dx * i;
        f1 = func(x1);
        integral += refine(func, f0, f1, x0, x1, err);
        x0 = x1;
        f0 = f1;
    }
    return integral;
}

double simpson_integral(double (*func)(double), double a, double b, int n)
{
    double dx = (b - a) / n;
    double fractor = dx / 6;
    double integral = 0;
    for (int i = 0; i < n; i++)
    {
        integral += fractor * (func(a + i * dx) + 4.0 * func(a + (0.5 + i) * dx) + func(a + (i + 1) * dx));
    }
    return integral;
}

double richardson_integral(double (*func)(double), double a, double b, int n, int k)
{
    double i_array[k + 1] = {0};
    int pow = 1;
    for (int i = 0; i < k + 1; i++)
    {
        i_array[i] = trapezoid_integral(func, a, b, n * pow);
        pow = pow * 2;
    }
    pow = 1;
    for (int i = 1; i < k + 1; i++)
    {
        pow = pow * 4;
        for (int j = 0; j < k - i + 1; j++)
        {
            i_array[j] = (i_array[j + 1] * pow - i_array[j]) / (pow - 1.0);
        }
    }
    return i_array[0];
}

double chebyshev_anomaly_integral(double (*func)(double), double a, double b, int n)
{
    double integral = 0;
    double x_i = 0;
    for (int i = 1; i < n + 1; i++)
    {
        x_i = (a + b) / 2 + (b - a) / 2 * cos(PI * (2 * i - 1) / (2 * n));
        integral += func(x_i);
    }
    return PI / n * integral;
}

double chebyshev_second_integral(double (*func)(double), double a, double b, int n)
{
    double integral = 0;
    double x_i = 0;
    for (int i = 1; i < n + 1; i++)
    {
        x_i = (a + b) / 2 + (b - a) / 2 * cos(PI * (2 * i - 1) / (2 * n));
        integral += func(x_i) * (1 - x_i * x_i);
    }
    return PI / n * integral;
}