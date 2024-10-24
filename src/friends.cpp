#include "s21_matrix_oop.h"

S21Matrix operator*(double num, const S21Matrix& other) {
  S21Matrix copy_matrix(other);
  copy_matrix.MulNumber(num);
  return copy_matrix;
}

S21Matrix operator*(const S21Matrix& other, double num) {
  S21Matrix copy_matrix(other);
  copy_matrix.MulNumber(num);
  return copy_matrix;
}
