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
    int N = 1024;  // the dimension of matrix
    My_Matrix H(N, N);
    My_Matrix R(N, N);
    My_Matrix Q(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = i - 1; j < i + 2; j++) {
        if (i == j) {
          H.set_element(i, j, 2);
        } else {
          H.set_element(i, (j + N) % N, -1);
        }
      }
    }
    long double v[N] = {0};
    int cnt = 0;
    cnt = H.Eigenvalue_sym(v, N, 1e-15);
    double u[N] = {0};
    for (int i = 0; i < N; i++){
      u[i] = (double)v[i];
    } ofstream out;
    out.open("1_eigenvalue_1024.dat", std::ofstream::binary);
    out.write(reinterpret_cast<const char *>(u), sizeof(double) * N);
    out.close();
    cout << "dimention of the matrix is " << N << endl;
    cout << "QR iteration steps is ";
    cout << cnt << endl;
    return 0;
  } catch (std::exception &e) {
    std::cerr << e.what() << endl;
  } catch (const char *error) {
    std::cerr << error << endl;
  }
  return 0;
}