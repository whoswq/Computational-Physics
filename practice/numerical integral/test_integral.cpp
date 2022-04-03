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

double quartic(double x)
{
    double x_2 = x * x;
    return x_2 * x_2;
}
double lorentz(double x)
{
    return 1.0 / (1.0 + x * x);
}
double exp_e(double x)
{
    return exp(-x);
}
double pendulum(double theta)
{
    return sqrt(2.0) * sqrt((PI * PI / 4 - theta * theta) / (cos(theta) - cos(PI / 2)));
}
// test second kind of chebyshev integral, f(x) = 1/4 x^4
double test_chebshev_2(double x){
    return quartic(x) * sqrt(1 - x * x);
}

int main()
{
    //cout << midpoint_integral(exp_e, 0.0, 5.0, 10) << endl;
    //cout << trapezoid_integral(exp_e, 0.0, 5.0, 10) << endl;
    //cout << self_adaptive_integral(exp_e, 0.0, 5.0, 10) << endl;
    //cout << simpson_integral(exp_e, 0.0, 5.0, 10) << endl;
    cout << richardson_integral(test_chebshev_2, -1, 1, 20, 2) << endl;
    cout << chebyshev_second_integral(quartic, -1, 1, 20) << endl;
    cout << chebyshev_anomaly_integral(pendulum, -PI / 2, PI / 2, 10) << endl;
    return 0;
}