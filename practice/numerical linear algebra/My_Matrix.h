#ifndef _My_Matrix_H_
#define _My_Matrix_H_

const double delta = 1e-15;

class My_Matrix {
 private:
  int row_number;
  int col_number;
  double **element;
  double determiniant;

 public:
  // constructor
  My_Matrix(int row, int col);
  // copying function
  My_Matrix(const My_Matrix &M);
  // destructor
  ~My_Matrix();
  // set one element of matrix
  void set_element(int row, int col, double ele);
  // set all element of matrix
  void set_matrix(double **ele);
  // return ij element of matrix
  double read_element(int i, int j) const;
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
  My_Matrix operator*(const double b);
  friend My_Matrix operator*(const double a, const My_Matrix &B);
  My_Matrix operator[](const My_Matrix &B);
  My_Matrix operator+(const My_Matrix &B);
  My_Matrix operator+(const double b);
  friend My_Matrix operator+(const double a, const My_Matrix &B);
  My_Matrix operator-(const My_Matrix &B);
  My_Matrix operator-(const double b);
  friend My_Matrix operator-(const double a, const My_Matrix &B);

  // calculate determinant
  double det_cal() const;
  double det_show() const;
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
                       const double eps) const;
  int Gauss_Seidel_iteration(My_Matrix &x, const My_Matrix &b, const int M,
                             const double eps) const;
  int Conjugate_Gradiant(My_Matrix &x, const My_Matrix &b,
                         const double eps) const;
  My_Matrix Gram_Schmidt_QR(My_Matrix &R);
};
// namespace My_Matrix
#endif