#ifndef S21_MATRIX_OOP_HPP
#define S21_MATRIX_OOP_HPP

#include <time.h>
#include <unistd.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include "eigen/Eigen/Dense"

#define MIN_DIM 1
#define MAX_DIM 46340
#define EPSILON 0.0001

class S21Matrix {
  friend S21Matrix operator*(double num, const S21Matrix& other);
  friend S21Matrix operator*(const S21Matrix& other, double num);

 private:
  int rows_, cols_;
  double* matrix_ = nullptr;

 public:
  int matrix_len = 0;

  /*
  //	Constructors & Destructors
  */
  S21Matrix();
  S21Matrix(int rows_, int cols_);
  S21Matrix(const S21Matrix& src);
  S21Matrix(S21Matrix&& src);
  ~S21Matrix();

  S21Matrix(Eigen::MatrixXd& src);  // for testing purpose

  /*
  //	Accessor & Mutator
  */
  int getRows();
  void editRows(int new_rows);
  int getCols();
  void editCols(int new_cols);
  double getValueByIndex(int index);
  double& setValueByIndex(int index);

  /*
  //	Functions-Crutches
  */
  void Describe();
  void RoundZeroes();
  int CustomPower(int init_value, int power);
  void CheckInverse();
  void AssignMatrix(const S21Matrix& other);
  void FillMatrix(int limit);
  void ShowMatrix();
  void Complements2x2(S21Matrix& minors);
  double Det2x2();
  S21Matrix CropMatrix(int i);

  /*
  //	Main Public Methods
  */
  bool EqMatrix(const S21Matrix& other);
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  double Determinant();
  S21Matrix Transpose();
  S21Matrix CalcComplements();
  S21Matrix Inverse();

  /*
  //	Operators Overloading
  */
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  bool operator==(const S21Matrix& other);
  void operator=(const S21Matrix& other);
  void operator+=(const S21Matrix& other);
  void operator-=(const S21Matrix& other);
  void operator*=(const S21Matrix& other);
  void operator*=(double num);
  double operator()(int i, int j);
};

#endif
