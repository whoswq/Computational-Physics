#ifndef _My_Matrix_H_
#define _My_Matrix_H_

const long double delta = 1e-15;
const long double PI = 3.1415926535897932384626433;

class My_Matrix {
 private:
  int row_number;
  int col_number;
  long double **element;
  long double determiniant;

 public:
  // constructor
  My_Matrix(int row, int col);
  // copying function
  My_Matrix(const My_Matrix &M);
  // destructor
  ~My_Matrix();
  // set one element of matrix
  void set_element(int row, int col, long double ele);
  // set all element of matrix
  void set_matrix(long double **ele);
  // return ij element of matrix
  long double read_element(int i, int j) const;
  // return row number
  int n_row() const;
  // return column number
  int n_col() const;
  // print the matrix
  void show_matrix() const;
  // joint two matrix in horizontal
  void joint_h(const My_Matrix &A, const My_Matrix &B);
  // joint two matrix in vertical
  void joint_v(const My_Matrix &A, const My_Matrix &B);
  int product(const My_Matrix &A, const My_Matrix &B);

  // reload operator
  My_Matrix &operator=(const My_Matrix &B);
  My_Matrix operator*(const My_Matrix &B);
  My_Matrix operator*(const long double b);
  friend My_Matrix operator*(const long double a, const My_Matrix &B);
  My_Matrix operator[](const My_Matrix &B);
  My_Matrix operator+(const My_Matrix &B);
  My_Matrix operator+(const long double b);
  friend My_Matrix operator+(const long double a, const My_Matrix &B);
  My_Matrix operator-(const My_Matrix &B);
  My_Matrix operator-(const long double b);
  friend My_Matrix operator-(const long double a, const My_Matrix &B);

  // calculate determinant
  long double det_cal() const;
  long double det_show() const;
  // interchange two rows
  void swap(int i, int j);
  // solve triangular matrix equation
  My_Matrix forward_L(const My_Matrix &B);
  My_Matrix backward_U(const My_Matrix &B);
  // store the results in U L
  My_Matrix LU_decomposition(My_Matrix &U, My_Matrix &L) const;
  // store the results in this
  My_Matrix LU_decomposition();
  My_Matrix LU_solvers(const My_Matrix &b) const;
  int Cholesky_decomposition(bool check_pos_def);
  int LDL_T_decomposition();
  My_Matrix Cholesky_decomposition(My_Matrix &L);
  My_Matrix Cholesky_solver(const My_Matrix &b) const;
  int Jacobi_iteration(My_Matrix &x, const My_Matrix &b, const int M,
                       const long double eps) const;
  int Gauss_Seidel_iteration(My_Matrix &x, const My_Matrix &b, const int M,
                             const long double eps) const;
  int Conjugate_Gradiant(My_Matrix &x, const My_Matrix &b,
                         const long double eps) const;
  int Gram_Schmidt_QR(My_Matrix &Q, My_Matrix &R) const;
  int Householder_QR(My_Matrix &Q, My_Matrix &R) const;
  int Givens_QR(My_Matrix &Q, My_Matrix &R) const;
  int Householder_Reduction(My_Matrix &H, My_Matrix &U) const;
  int Householder_Reduction_sym(My_Matrix &H) const;
  int upperHessenberg_QR(My_Matrix &Q, My_Matrix &R) const;
  int Triangular_QR(My_Matrix &Q, My_Matrix &R) const;
  int Wilkinson_QR_iter(int begin, int end);
  int Eigenvalue(long double *v, const int n) const;
  int Eigenvalue_sym(long double *v, const int n, const long double eps) const;
  int Eigenvalue_Triangular(long double *v, const int n, const long double eps) const;
};
int Householder(long double *x, long double *v, const int n, long double beta);
int Givens(long double alpha, long double beta, long double c, long double s, long double eta);
// namespace My_Matrix
#endif