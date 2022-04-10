#include "My_Matrix.h"

#include <math.h>
//#include <stdlib.h>

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
  this->element = new long double *[row];
  for (int i = 0; i < row; i++) {
    this->element[i] = new long double[col];
  }
  // set all elements to zero
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      this->element[i][j] = 0.0;
    }
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
  this->element = new long double *[M.n_row()];
  for (int i = 0; i < this->row_number; i++) {
    this->element[i] = new long double[this->col_number];
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

My_Matrix My_Matrix::operator+(const long double b) {
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

My_Matrix operator+(const long double a, const My_Matrix &B) {
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
      long double sum = 0;
      for (int k = 0; k < this->col_number; k++) {
        sum += this->element[i][k] * B.read_element(k, j);
      }
      C.set_element(i, j, sum);
    }
  }
  return C;
}

My_Matrix operator*(const long double a, const My_Matrix &B) {
  My_Matrix C(B.n_row(), B.n_col());
  for (int i = 0; i < B.n_row(); i++) {
    for (int j = 0; j < B.n_col(); j++) {
      C.set_element(i, j, B.read_element(i, j) * a);
    }
  }
  return C;
}

My_Matrix My_Matrix::operator*(const long double b) {
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

My_Matrix My_Matrix::operator-(const long double b) {
  return *this + (-1.0) * b;
}

My_Matrix operator-(const long double a, const My_Matrix &B) {
  return a + (-1.0) * B;
}

int My_Matrix::n_row() const { return row_number; }
int My_Matrix::n_col() const { return col_number; }

void My_Matrix::set_element(int row, int col, long double ele) {
  element[row][col] = ele;
}

void My_Matrix::set_matrix(long double **ele) {
  for (int i = 0; i < this->row_number; i++) {
    for (int j = 0; j < this->col_number; j++) {
      this->element[i][j] = ele[i][j];
    }
  }
}

long double My_Matrix::read_element(int i, int j) const {
  return element[i][j];
}

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
  long double *a = this->element[i];
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
      long double sum = 0;
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
    throw "My_Matrix, Error in joint_h, A && B have different row number";
  }
  if (this->col_number != (A.n_col() + B.n_col())) {
    throw "My_Matrix, Error in joint_h, this does not have enough columns";
  }
  for (int i = 0; i < A.n_row(); i++) {
    for (int j = 0; j < A.n_col(); j++) {
      this->element[i][j] = A.read_element(i, j);
    }
    for (int j = 0; j < B.n_col(); j++) {
      this->element[i][j + A.n_col()] = B.read_element(i, j);
    }
  }
}

void My_Matrix::joint_v(const My_Matrix &A, const My_Matrix &B) {
  if (A.n_col() != B.n_col()) {
    throw "My_Matrix, Error in joint_v, A && B have different column number";
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
  long double x = 0;
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
  long double x = 0;
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
  long double a_max = 0;
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

long double My_Matrix::det_cal() const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in det_cal, not a square matrix";
  }
  My_Matrix C = *this;
  C.LU_decomposition();
  long double det = 1;
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
  long double a_max = 0;
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
                                const long double eps = 1e-6) const {
  if (this->col_number != this->row_number) {
    throw "My_Matrix, Error in Jacobi_iteration, not a square Matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Jacobi_iteration, x b do not have same column number";
    return -1;
  }
  long double x_[this->col_number] = {0};
  long double norm = 0;
  for (int k = 0; k < M; k++) {
    if (k == M - 1) {
      throw "My_Matix, error in Jacobi_iteration, Not converged at given M && esp";
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
                                      const long double eps = 1e-6) const {
  if (this->col_number != this->row_number) {
    throw "My_Matrix, Error in Gauss_seidel_iteration, not a square matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Jacobi_iteration, x b do not have same column number";
    return -1;
  }
  long double norm = 0;
  for (int k = 0; k < M; k++) {
    if (k == M - 1) {
      throw "My_Matix, error in Jacobi_iteration, Not converged at given M && esp";
      return -1;
    }
    for (int col = 0; col < b.n_col(); col++) {
      long double x_new[this->row_number] = {0};
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
                                  const long double eps = 1e-6) const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Conjugate_Gradiant, not a square matrix";
    return -1;
  }
  if (x.n_col() != b.n_col()) {
    throw "My_Matrix, Error in Conjugate_Gradiant, not the same dimension";
    return -1;
  }
  // store residual vector
  long double r[this->row_number][x.n_col()] = {0};
  // store solution
  long double x_[this->row_number][x.n_col()] = {0};
  long double q[this->row_number][x.n_col()] = {0};
  long double alpha = 0;
  long double beta = 0;
  long double oldResi = 0;
  long double newResi = 0;
  long double pAp = 0;
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
  long double p[this->row_number][x.n_col()] = {0};
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
        long double var = 0;
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

int My_Matrix::Gram_Schmidt_QR(My_Matrix &Q, My_Matrix &R) const {
  if (Q.n_col() != this->col_number or Q.n_row() != this->row_number) {
    throw "My_Matrix, Error in Gram_Schmidt_QR, do not have proper dimension";
    return -1;
  }
  if (R.n_col() != R.n_row() or this->col_number != R.n_col()) {
    throw "My_Matrix, Error in Gram_Schmidt_QR, do not have proper dimension";
    return -1;
  }
  for (int i = 0; i < this->col_number; i++) {
    for (int j = 0; j < this->row_number; j++) {
      Q.set_element(j, i, this->read_element(j, i));
    }
    for (int j = 0; j < i; j++) {
      long double var = 0;
      for (int k = 0; k < this->row_number; k++) {
        var += Q.read_element(k, j) * Q.read_element(k, i);
      }
      R.set_element(j, i, var);
      for (int k = 0; k < this->row_number; k++) {
        var =
            Q.read_element(k, i) - R.read_element(j, i) * Q.read_element(k, j);
        Q.set_element(k, i, var);
      }
    }
    long double norm = 0;
    for (int k = 0; k < this->row_number; k++) {
      norm += Q.read_element(k, i) * Q.read_element(k, i);
    }
    R.set_element(i, i, sqrt(norm));
    if (norm < 1e-12) {
      return 0;
    }
    for (int k = 0; k < this->row_number; k++) {
      Q.set_element(k, i, Q.read_element(k, i) / R.read_element(i, i));
    }
  }
  return 0;
}

// this is not from textbook but from numerical linear algebra
// code in P122 of textbook may be not stable
int Householder(long double *x, long double *v, const int n,
                long double *beta) {
  long double eta = abs(x[0]);  // the maxium of x
  for (int i = 1; i < n; i++) {
    if (eta < abs(x[i])) {
      eta = abs(x[i]);
    }
  }
  for (int k = 0; k < n; k++) {
    x[k] = x[k] / eta;
  }
  long double sigma = 0;
  for (int k = 1; k < n; k++) {
    sigma += x[k] * x[k];
    v[k] = x[k];
  }
  if (sigma == 0) {
    beta[0] = 0;
  } else {
    long double alpha = sqrt(x[0] * x[0] + sigma);
    if (x[0] <= 0) {
      v[0] = x[0] - alpha;
    } else {
      v[0] = -sigma / (x[0] + alpha);
    }
    beta[0] = 2.0 * v[0] * v[0] / (sigma + v[0] * v[0]);
    long double v_0 = v[0];
    for (int k = 0; k < n; k++) {
      v[k] = v[k] / v_0;
    }
  }
  return 0;
}

int My_Matrix::Householder_QR(My_Matrix &Q, My_Matrix &R) const {
  if (Q.n_col() != this->col_number or Q.n_row() != this->row_number) {
    throw "My_Matrix, Error in Householder_QR, do not have proper dimension";
    return -1;
  }
  if (this->row_number != R.n_row() or this->col_number != R.n_row()) {
    throw "My_Matrxi, Error in Householder_QR, not the same dimension";
    return -1;
  }
  Q = *this;
  for (int j = 0; j < this->col_number; j++) {
    if (j < this->row_number) {
      long double v[this->row_number - j] = {0};
      long double beta[1] = {0};
      long double x[this->row_number - j] = {0};
      for (int k = j; k < this->row_number; k++) {
        x[k - j] = Q.read_element(k, j);
      }
      Householder(x, v, this->row_number - j, beta);
      long double vtA[this->row_number - j] = {0};
      for (int k = 0; k < this->col_number - j; k++) {
        long double var = 0;
        for (int row = 0; row < this->row_number - j; row++) {
          var += v[row] * Q.read_element(row + j, k + j);
        }
        vtA[k] = var;
      }
      for (int i = j; i < this->row_number; i++) {
        for (int k = j; k < this->col_number; k++) {
          long double var = 0;
          var = Q.read_element(i, k) - beta[0] * v[i - j] * vtA[k - j];
          Q.set_element(i, k, var);
        }
      }
      long double factor = sqrt(beta[0]);
      for (int k = 0; k < this->row_number; k++) {
        if (k <= j) {
          R.set_element(k, j, Q.read_element(k, j));
          Q.set_element(k, j, 0);
        } else {
          Q.set_element(k, j, v[k - j] * factor);
        }
      }
      Q.set_element(j, j, factor);
    }
  }
  return 0;
}

int Givens(long double alpha, long double beta, long double *c, long double *s,
           long double *eta) {
  eta[0] = sqrt(alpha * alpha + beta * beta);
  // if (eta[0] > 0) {
  //   c[0] = alpha / eta[0];
  //   s[0] = beta / eta[0];
  // } else {
  //   c[0] = 1;
  //   s[0] = 0;
  // }
  if (abs(beta) == 0) {
    c[0] = 1.0;
    s[0] = 0.0;
  } else {
    if (abs(beta) > abs(alpha)) {
      long double tau = alpha / beta;
      s[0] = 1.0 / sqrt(1.0 + tau * tau);
      c[0] = s[0] * tau;
    } else {
      long double tau = beta / alpha;
      c[0] = 1.0 / sqrt(1.0 + tau * tau);
      s[0] = c[0] * tau;
    }
  }
  return 0;
}

/*
store c, s in Q
*/
int My_Matrix::Givens_QR(My_Matrix &Q, My_Matrix &R) const {
  if (Q.n_col() != this->col_number or Q.n_row() != this->row_number) {
    throw "My_Matrix, Error in Givens_QR, do not have proper dimension";
    return -1;
  }
  if (this->row_number != R.n_row() or this->col_number != R.n_row()) {
    throw "My_Matrxi, Error in Givens_QR, not the same dimension";
    return -1;
  }
  R = *this;
  long double c[1] = {0};
  long double s[1] = {0};
  long double eta[1] = {0};
  for (int i = 0; i < this->row_number - 1; i++) {
    for (int j = i + 1; j < this->row_number; j++) {
      Givens(R.read_element(i, i), R.read_element(j, i), c, s, eta);
      Q.set_element(i, j, c[0]);
      Q.set_element(j, i, s[0]);
      if (eta[0] > 1e-50) {
        R.set_element(i, i, eta[0]);
        R.set_element(j, i, 0);
        for (int k = i + 1; k < R.n_col(); k++) {
          long double u = R.read_element(i, k);
          long double v = R.read_element(j, k);
          R.set_element(i, k, c[0] * u + s[0] * v);
          R.set_element(j, k, -s[0] * u + c[0] * v);
        }
      }
    }
  }
  return 0;
}

int My_Matrix::Householder_Reduction(My_Matrix &H, My_Matrix &U) const {
  return 0;
  // to do
}

/*
reduce a real symmetric matrix to upper Hessenberg Matrix
*/
int My_Matrix::Householder_Reduction_sym(My_Matrix &H) const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Householder_Reduction_sym, not a square matrix";
    return -1;
  }
  H = *this;
  for (int k = 0; k < this->col_number - 2; k++) {
    long double v[this->row_number - k - 1] = {0};
    long double beta[1] = {0};
    long double x[this->row_number - k - 1] = {0};
    for (int j = k + 1; j < this->row_number; j++) {
      x[j - k - 1] = H.read_element(j, k);
    }
    Householder(x, v, this->row_number - k - 1, beta);
    long double u[this->row_number - k - 1] = {0};
    for (int i = k + 1; i < this->row_number; i++) {
      long double var = 0;
      for (int j = k + 1; j < this->row_number; j++) {
        var += H.read_element(i, j) * v[j - k - 1];
      }
      u[i - k - 1] = var * beta[0];
    }
    long double w[this->row_number - k - 1] = {0};
    long double var = 0;
    for (int i = 0; i < this->row_number - k - 1; i++) {
      var += u[i] * v[i];
    }
    var = var * beta[0] / 2.0;
    for (int i = 0; i < this->row_number - k - 1; i++) {
      w[i] = u[i] - var * v[i];
    }
    long double norm = 0;
    for (int i = k + 1; i < this->row_number; i++) {
      norm += H.read_element(i, k) * H.read_element(i, k);
    }
    H.set_element(k + 1, k, sqrt(norm));
    H.set_element(k, k + 1, H.read_element(k + 1, k));
    for (int i = k + 1; i < this->row_number; i++) {
      for (int j = k + 1; j < this->row_number; j++) {
        long double var = H.read_element(i, j) - v[i - k - 1] * w[j - k - 1] -
                          w[i - k - 1] * v[j - k - 1];
        H.set_element(i, j, var);
      }
    }
    // other off-diagonal term is 0
    for (int i = k + 2; i < this->row_number; i++) {
      H.set_element(i, k, 0.0);
      H.set_element(k, i, 0.0);
    }
  }
  return 0;
}

/*

*/
int My_Matrix::upperHessenberg_QR(My_Matrix &Q, My_Matrix &R) const {
  return 0;
  // to do
}

/*
QR decomposition of a triangular matrix
*/
int My_Matrix::Triangular_QR(My_Matrix &Q, My_Matrix &R) const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Triangular_QR, not a square matrix";
    return -1;
  }
  R = *this;
  long double c[1] = {0};
  long double s[1] = {0};
  long double eta[1] = {0};
  for (int i = 0; i < this->row_number - 1; i++) {
    Givens(R.read_element(i, i), R.read_element(i + 1, i), c, s, eta);
    Q.set_element(i, i + 1, c[0]);
    Q.set_element(i + 1, i, s[0]);
    if (eta[0] > 1e-50) {
      R.set_element(i, i, eta[0]);
      R.set_element(i + 1, i, 0);
      for (int k = i + 1; k < i + 3; k++) {
        if (k < R.n_col()) {
          long double u = R.read_element(i, k);
          long double v = R.read_element(i + 1, k);
          R.set_element(i, k, c[0] * u + s[0] * v);
          R.set_element(i + 1, k, -s[0] * u + c[0] * v);
        }
      }
    }
  }
  return 0;
}

/*
need not to calculate rotation matrix
*/
int My_Matrix::Wilkinson_QR_iter(int start, int final) {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Wilkinson_QR_iter, not a square matrix";
    return -1;
  }
  if (start >= this->row_number or final >= this->row_number) {
    throw "My_Matrix, Error in Wilkinson_QR_iter, index out of range";
    return -1;
  }
  long double d = (this->read_element(final - 1, final - 1) -
                   this->read_element(final, final)) /
                  2.0;
  long double u = this->read_element(final, final - 1) *
                  this->read_element(final, final - 1);
  long double mu = 0;
  if (abs(d) == 0) {
    // may be matrix has two same diagonal terms
    mu = this->read_element(final, final) + sqrt(u);
  } else {
    mu = this->read_element(final, final) -
         u / (d + d / abs(d) * sqrt(d * d + u));
  }
  long double x = this->read_element(start, start) - mu;
  long double z = this->read_element(start + 1, start);
  long double c[1] = {0};
  long double s[1] = {0};
  long double eta[1] = {0};
  for (int k = start; k < final; k++) {
    Givens(x, z, c, s, eta);
    // calculate G@T; row: k, k+1; col: k-1, k, k+1, k+2
    for (int i = k - 1; i < k + 3; i++) {
      if (i >= start && i <= final) {
        long double u = this->read_element(k, i);
        long double v = this->read_element(k + 1, i);
        this->set_element(k, i, c[0] * u + s[0] * v);
        this->set_element(k + 1, i, -s[0] * u + c[0] * v);
      }
    }

    // calculate G@T@G^t; row: k-1, k, k+1, k+2; col: k, k+1
    for (int i = k - 1; i < k + 3; i++) {
      if (i >= start && i <= final) {
        long double u = this->read_element(i, k);
        long double v = this->read_element(i, k + 1);
        this->set_element(i, k, u * c[0] + v * s[0]);
        this->set_element(i, k + 1, -u * s[0] + v * c[0]);
      }
    }
    if (k <= final - 2) {
      x = this->read_element(k + 1, k);
      z = this->read_element(k + 2, k);
    }
  }
  return 0;
}

/*
Eigen value of a genenral matrix
*/
int My_Matrix::Eigenvalue(long double *v, const int n) const {
  return 0;
  // todo
}

/*
Eigen value of a symmetric matrix
*/
int My_Matrix::Eigenvalue_sym(long double *v, const int n,
                              const long double eps = 1e-15) const {
  if (this->row_number != this->col_number) {
    throw "My_Matrix, Error in Eigenvalue_sym, not a square matrix";
    return -1;
  }
  My_Matrix H(this->n_row(), this->n_col());
  this->Householder_Reduction_sym(H);  // reduce to trangular matrix
  // if n_col() is large, this step uses much more time than iteration
  cout << "finished Householder reduction" << endl;
  int final = H.n_col() - 1;
  int cnt = 0;
  while (true) {
    long double f = final;
    // check whether H is digonal matrix;
    for (int k = f; k > 0; k--) {
      if (abs(H.read_element(k, k - 1)) <
          eps *
              (abs(H.read_element(k, k)) + abs(H.read_element(k - 1, k - 1)))) {
        H.set_element(k, k - 1, 0);
        H.set_element(k - 1, k, 0);
        final = k - 1;
      } else {
        final = k;
        break;
      }
    }
    // one step of QR iteration
    if (final == 0) {
      break;
    }
    // cout << "final = " << final << endl;
    // cout << "cnt = " << cnt << endl;
    H.Wilkinson_QR_iter(0, final);
    // cout << cnt << endl;
    cnt += 1;
  }
  for (int k = 0; k < n; k++) {
    v[k] = H.read_element(k, k);
  }
  return cnt;
}

/*
Eigen value of triangular matrix
*/
int My_Matrix::Eigenvalue_Triangular(long double *v, const int n,
                                     const long double eps = 1e-15) const {
  My_Matrix H = *this;
  int start = 0;
  int final = H.n_col() - 1;
  while (true) {
    long double s = start;
    long double f = final;
    // check whether H is digonal matrix;
    for (int k = s; k < f; k++) {
      if (abs(H.read_element(k, k + 1)) <
          eps *
              (abs(H.read_element(k, k)) + abs(H.read_element(k + 1, k + 1)))) {
        H.set_element(k, k + 1, 0);
        H.set_element(k + 1, k, 0);
      } else {
        start = k;
        break;
      }
    }
    for (int k = f; k > s; k--) {
      if (abs(H.read_element(k, k - 1)) <
          eps *
              (abs(H.read_element(k, k)) + abs(H.read_element(k - 1, k - 1)))) {
        H.set_element(k, k - 1, 0);
        H.set_element(k - 1, k, 0);
        final = k - 1;
      } else {
        final = k;
        break;
      }
    }
    // one step of QR iteration
    if (start == final) {
      break;
    }
    H.Wilkinson_QR_iter(start, final);
  }
  for (int k = 0; k < n; k++) {
    v[k] = H.read_element(k, k);
  }
  return 0;
}

/*
int main() {
  try {
    int N = 2000;  // 370, 371
    long double gama = 2;
    long double alpha = 0.02;
    My_Matrix H(N, N);
    My_Matrix R(N, N);
    My_Matrix Q(N, N);
    long double a[9] = {0, 1, 1, 0};
    for (int i = 0; i < N; i++) {
      for (int j = i - 1; j < i + 2; j++) {
        if (i == j) {
          H.set_element(i, j, 2.0 + gama * cos(2.0 * PI * j * alpha));
        } else {
          H.set_element(i, (j + N) % N, -1);
        }
      }
    }
    long double v[N] = {0};
    int cnt = 0;
    cnt = H.Eigenvalue_sym(v, N, 1e-12);
    double u[N] = {0};
    for (int i = 0; i < N; i++) {
      u[i] = (double)v[i];
    }
    cout << "total iteration steps is ";
    cout << cnt << endl;
    ofstream out;
    out.open("test.dat", std::ofstream::binary);
    out.write(reinterpret_cast<const char *>(u), sizeof(double) * N);
    out.close();
    return 0;
  } catch (std::exception &e) {
    std::cerr << e.what() << endl;
  } catch (const char *error) {
    std::cerr << error << endl;
  }
}
*/