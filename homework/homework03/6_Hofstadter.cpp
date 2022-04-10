#include <math.h>
#include <stdlib.h>

#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "My_Matrix.h"
// using namespace std;

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

int main() {
  try {
    int N = 1000;  // the dimension of matrix
    int M = 1000;  // the number of alpha
    long double gama = 2;
    long double alpha = 0;
    long double da = (long double)1.0 / M;
    ofstream out;
    out.open("6_Hofstadter_1000_1000.dat", std::ofstream::binary);
    for (int m = 0; m < M; m++) {
      My_Matrix H(N, N);
      for (int i = 0; i < N; i++) {
        for (int j = i - 1; j < i + 2; j++) {
          if (i == j) {
            H.set_element(i, j, 2.0 + gama * cos(2.0 * PI * j * alpha));
          } else {
            H.set_element(i, (j + N) % N, -1);
          }
        }
      }
      alpha += da;
      long double v[N] = {0};
      int cnt = 0;
      cnt = H.Eigenvalue_sym(v, N, 1e-15);
      double u[N] = {0};
      for (int i = 0; i < N; i++) {
        u[i] = (double)v[i];
      }
      out.write(reinterpret_cast<const char *>(u), sizeof(double) * N);
      cout << "dimension of the matrix is " << N << endl;
      cout << "QR iteration steps is ";
      cout << cnt << endl;
      cout << "total matrix number is " << M << ", "
           << " now is " << m << endl;
    }
    out.close();
    return 0;
  } catch (std::exception &e) {
    std::cerr << e.what() << endl;
  } catch (const char *error) {
    std::cerr << error << endl;
  }
  return 0;
}