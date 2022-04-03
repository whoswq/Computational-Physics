#include "My_Matrix.h"

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
// using namespace std;

using std::cin;
using std::cout;
using std::endl;
using std::ofstream;
using std::string;
using std::stringstream;
using std::to_string;

My_Matrix::My_Matrix(int row, int col) {
  this->row_number = row;
  this->col_number = col;
  this->element = new double *[row];
  for (int i = 0; i < row; i++) {
    this->element[i] = new double[col];
  }
}

My_Matrix::~My_Matrix() {
  for (int i = 0; i < row_number; i++) {
    // delete[] this->element[i];
  }
  // delete[] this->element;
  this->element = NULL;
  this->row_number = 0;
  this->col_number = 0;
}

My_Matrix::My_Matrix(const My_Matrix &M) {
  this->row_number = M.n_row();
  this->col_number = M.n_col();
  this->element = new double *[M.n_row()];
  for (int i = 0; i < this->row_number; i++) {
    this->element[i] = new double[this->col_number];
  }
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      this->element[i][j] = M.read_element(i, j);
    }
  }
}

My_Matrix &My_Matrix::operator=(const My_Matrix &B) {
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      this->set_element(i, j, B.read_element(i, j));
    }
  }
  return *this;
}

My_Matrix My_Matrix::operator+(const My_Matrix &B) {
  if (this->row_number != B.n_row() or this->col_number != B.n_col()) {
    throw "My_Matrix, Error in operator+";
  }
  My_Matrix C(this->row_number, this->col_number);
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      C.set_element(i, j, this->read_element(i, j) + B.read_element(i, j));
    }
  }
  return C;
}

My_Matrix My_Matrix::operator+(const double b) {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in operator+";
  }
  My_Matrix C(this->row_number, this->col_number);
  C = *this;
  for (int k = 0; k < this->row_number; k++) {
    C.set_element(k, k, this->element[k][k] + b);
  }
  return C;
}

My_Matrix operator+(const double a, const My_Matrix &B) {
  if (B.n_row() != B.n_row()) {
    throw "My_Matrix, Error in operator+";
  }
  My_Matrix C(B.n_row(), B.n_col());
  C = B;
  for (int k = 0; k < B.n_row(); k++) {
    C.set_element(k, k, B.read_element(k, k) + a);
  }
  return C;
}

My_Matrix My_Matrix::operator*(const My_Matrix &B) {
  if (this->col_number != B.n_row()) {
    throw "My_Matrix, Error in operator*";
  }
  My_Matrix C(this->row_number, B.n_col());
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < B.n_col(); j++) {
      double sum = 0;
      for (int k = 0; k < this->col_number; k++) {
        sum += this->element[i][k] * B.read_element(k, j);
      }
      C.set_element(i, j, sum);
    }
  }
  return C;
}

My_Matrix operator*(const double a, const My_Matrix &B) {
  My_Matrix C(B.n_row(), B.n_col());
  for (int i = 0; i < B.n_row(); i++) {
    for (int j = 0; j < B.n_col(); j++) {
      C.set_element(i, j, B.read_element(i, j) * a);
    }
  }
  return C;
}

My_Matrix My_Matrix::operator*(const double b) {
  My_Matrix C(this->row_number, this->col_number);
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      C.set_element(i, j, this->element[i][j] * b);
    }
  }
  return C;
}

My_Matrix My_Matrix::operator-(const My_Matrix &B) {
  return *this + (-1.0) * B;
}

My_Matrix My_Matrix::operator-(const double b) { return *this + (-1.0) * b; }

My_Matrix operator-(const double a, const My_Matrix &B) {
  return a + (-1.0) * B;
}

int My_Matrix::n_row() const { return row_number; }
int My_Matrix::n_col() const { return col_number; }

void My_Matrix::set_element(int row, int col, double ele) {
  element[row][col] = ele;
}

void My_Matrix::set_matrix(double **ele) {
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      this->element[i][j] = ele[i][j];
    }
  }
}

double My_Matrix::read_element(int i, int j) const { return element[i][j]; }

void My_Matrix::show_matrix() const {
  cout << "[" << endl;
  for (int i = 0; i < row_number; i++) {
    for (int j = 0; j < col_number; j++) {
      cout << std::setiosflags(std::ios::scientific | std::ios::showpos)
           << element[i][j] << " ";
    }
    cout << endl;
  }
  cout << " ]" << endl;
}

void My_Matrix::swap(int i, int j) {
  if (i > this->row_number or j > this->row_number) {
    throw "My_Matrxi, Error in swap";
  }
  double *a = this->element[i];
  this->element[i] = this->element[j];
  this->element[j] = a;
}

int My_Matrix::product(const My_Matrix &A, const My_Matrix &B) {
  if (A.n_col() != B.n_row() or this->row_number != A.n_row() or
      this->col_number != B.n_col()) {
    throw "My_Matrix, error in product";
  }
  for (int i = 0; i < A.n_row(); i++) {
    for (int j = 0; j < B.n_col(); j++) {
      double sum = 0;
      for (int k = 0; k < A.n_col(); k++) {
        sum += A.read_element(i, k) * B.read_element(k, j);
      }
      this->element[i][j] = sum;
    }
  }
  return 0;
}

void My_Matrix::joint_h(const My_Matrix &A, const My_Matrix &B) {
  if (A.n_row() != B.n_row()) {
    throw "My_Matrix, Error in joint_h, A and B have different row number";
  }
  if (this->col_number != (A.n_col() + B.n_col())) {
    throw "My_Matrix, Error in joint_h, this does not have enough columns";
  }
  for (int i = 0; i < A.n_row(); i++) {
    for (int j = 0; j < A.n_col(); j++) {
      this->element[i][j] = A.read_element(i, j);
      this->show_matrix();
    }
    for (int j = 0; j < B.n_col(); j++) {
      this->element[i][j + A.n_col()] = B.read_element(i, j);
    }
  }
}

void My_Matrix::joint_v(const My_Matrix &A, const My_Matrix &B) {
  if (A.n_col() != B.n_col()) {
    throw "My_Matrix, Error in joint_v, A and B have different column number";
  }
  if (this->row_number != (A.n_row() + B.n_row())) {
    throw "My_Matrix, Error in joint_v, this does not have enough rows";
  }
  for (int j = 0; j < A.n_col(); j++) {
    for (int i = 0; i < A.n_row(); i++) {
      this->element[i][j] = A.read_element(i, j);
    }
    for (int i = 0; i < B.n_row(); i++)
      this->element[i + A.n_row()][j] = B.read_element(i, j);
  }
}

My_Matrix My_Matrix::forward_L(const My_Matrix &B) {
  My_Matrix C(B.n_row(), B.n_col());
  double x = 0;
  for (int k = 0; k < B.n_col(); k++) {
    for (int i = 0; i < this->row_number; i++) {
      if (abs(this->element[i][i]) < delta) {
        throw "My_Matrix, Error in forward_L, this does not have inverse";
      }
      x = B.read_element(i, k);
      for (int j = 0; j < i; j++) {
        x -= this->element[i][j] * C.read_element(j, k);
      }
      C.set_element(i, k, x / this->element[i][i]);
    }
  }
  return C;
}

My_Matrix My_Matrix::backward_U(const My_Matrix &B) {
  My_Matrix C(B.n_row(), B.n_col());
  double x = 0;
  for (int k = 0; k < B.n_col(); k++) {
    for (int i = this->row_number - 1; i >= 0; i--) {
      if (abs(this->element[i][i]) < delta) {
        throw "My_Matrix, Error in backward_U, this does not have inverse";
      }
      x = B.read_element(i, k);
      for (int j = this->row_number - 1; j > i; j--) {
        x -= this->element[i][j] * C.read_element(j, k);
      }
      C.set_element(i, k, x / this->element[i][i]);
    }
  }
  return C;
}

My_Matrix My_Matrix::LU_decomposition() {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in LU_decomposition, not a square matrix";
  }
  double a_max = 0;
  int s = 0;
  My_Matrix P(this->row_number, 1);  // record the permutation
  for (int i = 0; i < this->row_number; i++) {
    P.set_element(i, 0, i);
  }
  for (int k = 0; k < this->row_number - 1; k++) {
    a_max = this->element[k][k];
    s = k;
    for (int i = k + 1; i < this->row_number; i++) {
      if (abs(a_max) < abs(this->element[i][k])) {
        a_max = this->element[i][k];
        s = i;  // select the max column element
      }
    }
    if (s != k) {
      this->swap(k, s);
      P.swap(k, s);
    }
    if (abs(this->element[k][k]) < delta) {
      throw "My_Matrix, Error in LU_decoposition, singular matrix";
    }
    for (int i = k + 1; i < this->row_number; i++) {
      this->element[i][k] = this->element[i][k] / this->element[k][k];  // l_ik
      for (int j = k + 1; j < this->row_number; j++) {
        this->element[i][j] =
            this->element[i][j] - this->element[i][k] * this->element[k][j];
      }
    }
  }
  return P;
}

double My_Matrix::det_cal() const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in det_cal, not a square matrix";
  }
  My_Matrix C = *this;
  C.LU_decomposition();
  double det = 1;
  for (int i = 0; i < this->col_number; i++) {
    det = det * C.read_element(i, i);
  }
  return det;
}

My_Matrix My_Matrix::LU_decomposition(My_Matrix &U, My_Matrix &L) const {
  if (this->col_number != U.n_col() or this->row_number != U.n_row()) {
    throw "My_Matrix, Error in LU_decomposition, not same size";
  }
  U = *this;
  My_Matrix P(this->row_number, 1);  // record the permutation
  for (int i = 0; i < this->row_number; i++) {
    P.set_element(i, 0, i);
  }
  double a_max = 0;
  int s = 0;
  for (int k = 0; k < U.n_row() - 1; k++) {
    a_max = U.read_element(k, k);
    s = k;
    for (int i = k + 1; i < this->row_number; i++) {
      if (abs(a_max) < abs(U.read_element(i, k))) {
        a_max = U.read_element(i, k);
        s = i;  // select the max column element
      }
    }
    if (s != k) {
      U.swap(k, s);
      L.swap(k, s);
      P.swap(k, s);
    }
    if (U.read_element(k, k) < delta) {
      throw "My_Matrix, Error in LU_decoposition, singular matrix";
    }
    L.set_element(k, k, 1.0);
    for (int i = k + 1; i < this->row_number; i++) {
      L.set_element(i, k, U.read_element(i, k) / U.read_element(k, k));  // l_ik
      for (int j = k + 1; j < this->row_number; j++) {
        U.set_element(
            i, j,
            U.read_element(i, j) - L.read_element(i, k) * U.read_element(k, j));
      }
    }
  }
  L.set_element(this->row_number - 1, this->row_number - 1, 1.0);
  for (int i = 0; i < this->row_number; i++) {
    for (int k = i + 1; k < this->row_number; k++) {
      L.set_element(i, k, 0);
      U.set_element(k, i, 0);
      cout << i << ' ' << k << endl;
    }
  }
  return P;
}

My_Matrix My_Matrix::LU_solvers(const My_Matrix &b) const {
  My_Matrix U(this->row_number, this->col_number);
  My_Matrix L(this->row_number, this->col_number);
  My_Matrix P(this->row_number, 1);
  My_Matrix PB(this->row_number, b.n_col());
  P = this->LU_decomposition(U, L);
  for (int k = 0; k < b.n_col(); k++) {
    for (int i = 0; i < b.n_row(); i++) {
      PB.set_element(i, k, b.read_element(P.read_element(i, 0), k));
    }
  }
  My_Matrix y = L.forward_L(PB);
  My_Matrix x = U.backward_U(y);
  return x;
}

int My_Matrix::Cholesky_decomposition(bool check_pos_def) {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Cholesky_decomposition, not a square matrix";
  }
  for (int k = 0; k < this->row_number; k++) {
    if (this->element[k][k] > 0) {
      this->element[k][k] = sqrt(this->element[k][k]);
    } else {
      check_pos_def = false;
      throw "My_Matrix, Error in Cholesky_decoposition, not positive difinite";
    }
    for (int i = k + 1; i < this->row_number; i++) {
      this->element[i][k] = this->element[i][k] / this->element[k][k];
    }
    for (int j = k + 1; j < this->row_number; j++) {
      for (int r = j; r < this->row_number; r++) {
        this->element[r][j] =
            this->element[r][j] - this->element[r][k] * this->element[j][k];
      }
    }
  }
  return 0;
}

My_Matrix My_Matrix::Cholesky_solver(const My_Matrix &b) const {
  My_Matrix B = *this;
  bool flag = true;
  B.Cholesky_decomposition(flag);
  for (int i = 0; i < this->row_number; i++) {
    for (int j = i + 1; j < this->row_number; j++) {
      B.set_element(i, j, B.read_element(j, i));
    }
  }
  My_Matrix y = B.forward_L(b);
  My_Matrix x = B.backward_U(y);
  return x;
}

int My_Matrix::Jacobi_iteration(My_Matrix &x, const My_Matrix &b,
                                const int M = 1e5,
                                const double eps = 1e-6) const {
  if (this->col_number != this->row_number) {
    throw "My_Matrix, Error in Jacobi_iteration, not a square Matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Jacobi_iteration, x b do not have same column number";
    return -1;
  }
  double x_[this->col_number] = {0};
  double norm = 0;
  for (int k = 0; k < M; k++) {
    if (k == M - 1) {
      throw "My_Matix, error in Jacobi_iteration, Not converged at given M and esp";
      return -1;
    }
    for (int col = 0; col < b.n_col(); col++) {
      for (int i = 0; i < this->col_number; i++) {
        x_[i] = b.read_element(i, col);
        for (int j = 0; j < this->col_number; j++) {
          if (i != j) {
            x_[i] = x_[i] - this->read_element(i, j) * x.read_element(j, col);
          }
        }
        x_[i] = x_[i] / this->read_element(i, i);
      }
      for (int k = 0; k < this->col_number; k++) {
        norm +=
            (x_[k] - x.read_element(k, col)) * (x_[k] - x.read_element(k, col));
      }
      if (sqrt(norm) < eps) {
        return 0;
      }
      for (int i = 0; i < this->col_number; i++) {
        x.set_element(i, col, x_[i]);
      }
      norm = 0;
    }
  }
  return 0;
}

int My_Matrix::Gauss_Seidel_iteration(My_Matrix &x, const My_Matrix &b,
                                      const int M = 1e5,
                                      const double eps = 1e-6) const {
  if (this->col_number != this->row_number) {
    throw "My_Matrix, Error in Gauss_seidel_iteration, not a square matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Jacobi_iteration, x b do not have same column number";
    return -1;
  }
  double norm = 0;
  for (int k = 0; k < M; k++) {
    if (k == M - 1) {
      throw "My_Matix, error in Jacobi_iteration, Not converged at given M and esp";
      return -1;
    }
    for (int col = 0; col < b.n_col(); col++) {
      double x_new[this->row_number] = {0};
      for (int row = 0; row < this->row_number; row++) {
        x_new[row] = x.read_element(row, col);
      }
      for (int i = 0; i < this->col_number; i++) {
        x_new[i] = b.read_element(i, col);
        for (int j = 0; j < this->row_number; j++) {
          if (i != j) {
            x_new[i] = x_new[i] - this->read_element(i, j) * x_new[j];
          }
        }
        x_new[i] = x_new[i] / this->read_element(i, i);
      }
      for (int i = 0; i < this->row_number; i++) {
        norm = (x_new[i] - x.read_element(i, col)) *
               (x_new[i] - x.read_element(i, col));
      }
      if (sqrt(norm) < eps) {
        return 0;
      }
      norm = 0;
      for (int i = 0; i < this->row_number; i++) {
        x.set_element(i, col, x_new[i]);
      }
    }
  }
  return 0;
}

int My_Matrix::Conjugate_Gradiant(My_Matrix &x, const My_Matrix &b,
                                  const double eps = 1e-6) const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Conjugate_Gradiant, not a square matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Conjugate_Gradiant, not the same dimension";
    return -1;
  }
  // store residual vector
  double r[this->row_number][x.n_col()] = {0};
  // store solution
  double x_[this->row_number][x.n_col()] = {0};
  double q[this->row_number][x.n_col()] = {0};
  double alpha = 0;
  double beta = 0;
  double oldResi = 0;
  double newResi = 0;
  double pAp = 0;
  for (int col = 0; col < x.n_col(); col++) {
    for (int i = 0; i < this->row_number; i++) {
      x_[i][col] = x.read_element(i, col);
    }
  }
  for (int col = 0; col < x.n_col(); col++) {
    for (int i = 0; i < this->row_number; i++) {
      r[i][col] = b.read_element(i, col);
      for (int j = 0; j < this->row_number; j++) {
        r[i][col] -= this->read_element(i, j) * x.read_element(j, col);
      }
    }
  }
  double p[this->row_number][x.n_col()] = {0};
  for (int col = 0; col < x.n_col(); col++) {
    for (int i = 0; i < this->row_number; i++) {
      p[i][col] = r[i][col];
    }
  }
  for (int col = 0; col < x.n_col(); col++) {
    oldResi = 0;
    for (int i = 0; i < this->col_number; i++) {
      oldResi += r[i][col] * r[i][col];
    }
    while (sqrt(oldResi) > eps) {
      for (int i = 0; i < this->row_number; i++) {
        double var = 0;
        for (int j = 0; j < this->row_number; j++) {
          var += this->read_element(i, j) * p[j][col];
        }
        q[i][col] = var;
      }
      for (int i = 0; i < this->row_number; i++) {
        pAp += p[i][col] * q[i][col];
      }
      alpha = oldResi / pAp;
      pAp = 0;
      for (int i = 0; i < this->row_number; i++) {
        x_[i][col] = x_[i][col] + alpha * p[i][col];
      }
      for (int i = 0; i < this->row_number; i++) {
        r[i][col] = r[i][col] - alpha * q[i][col];
      }
      newResi = 0;
      for (int i = 0; i < this->row_number; i++) {
        newResi += r[i][col] * r[i][col];
      }
      beta = newResi / oldResi;
      oldResi = newResi;
      for (int i = 0; i < this->row_number; i++) {
        p[i][col] = r[i][col] + beta * p[i][col];
      }
    }
  }
  for (int col = 0; col < x.n_col(); col++) {
    for (int i = 0; i < this->row_number; i++) {
      x.set_element(i, col, x_[i][col]);
    }
  }
  return 0;
}

My_Matrix My_Matrix::Gram_Schmidt_QR(My_Matrix &R) {
  My_Matrix Q(this->row_number, this->col_number);
  if (R.n_col() != R.n_row() or this->col_number != R.n_col()) {
    throw "My_Matrix, Error in Gram_Schmidt_QR, do not have proper dimension";
    return -1;
  }
  for (int i = 0; i < this->col_number; i++) {
    for (int j = 0; j < this->row_number; j++) {
      Q.set_element(j, i, this->read_element(j, i));
    }
    for (int j = 0; j < i - 1) {
      double var = 0;
      for (int k = 0; k < this->row_number; k++) {
        var += Q.read_element(k, j) * Q.read_element(k.i);
      }
      R.set_element(j, i, var);
      for (int k = 0; k < this->row_number; k++) {
        var =
            Q.read_element(k, i) - R.read_element(j, i) * Q.read_element(k, j);
        Q.set_element(k, i, var);
      }
    }
    double norm = 0;
    for (int k = 0; k < this->row_number; k++) {
      norm += Q.read_element(k, i) * Q.read_element(k, i);
    }
    R.set_element(i, i, sqrt(norm));
  }
}
int main() {
  try {
    My_Matrix A(4, 4);
    My_Matrix B(4, 1);
    My_Matrix x(4, 1);
    B.set_element(0, 0, 0);
    B.set_element(1, 0, 0);
    B.set_element(2, 0, 1);
    B.set_element(3, 0, 1);
    x.set_element(0, 0, 0);
    x.set_element(1, 0, 0);
    x.set_element(2, 0, 0);
    x.set_element(3, 0, 0);
    double a[16] = {4, -1, -1, 0, -1, 4, 0, -1, -1, 0, 4, -1, 0, -1, -1, 4};
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        A.set_element(i, j, a[i * 4 + j]);
      }
    }
    x.show_matrix();
    A.Conjugate_Gradiant(x, B);
    x.show_matrix();
    return 0;
  } catch (std::exception &e) {
    std::cerr << e.what() << endl;
  } catch (const char *error) {
    std::cerr << error << endl;
  }
}