#include "s21_matrix_oop.h"

S21Matrix::S21Matrix() {
  //	std::cout << "Constractor: [" << this << "]" << std::endl;
  rows_ = 2;
  cols_ = 2;
  matrix_len = rows_ * cols_;
  matrix_ = new double[matrix_len]();
}

S21Matrix::S21Matrix(int rows_, int cols_) : rows_(rows_), cols_(cols_) {
  if (rows_ < MIN_DIM || rows_ > MAX_DIM || cols_ < MIN_DIM ||
      cols_ > MAX_DIM) {
    throw std::invalid_argument(
        "Error in parameterized constructor: ∀rows_,cols_ ∈ [1, 46340]\nExit "
        "Failure\n");
    exit(1);
  }
  //	std::cout << "Constructor: [" << this << "]" << std::endl;
  matrix_len = rows_ * cols_;
  matrix_ = new double[matrix_len]();
}

S21Matrix::S21Matrix(const S21Matrix& src)
    : rows_(src.rows_), cols_(src.cols_) {
  matrix_len = rows_ * cols_;
  matrix_ = new double[matrix_len]();
  std::memcpy(matrix_, src.matrix_, matrix_len * sizeof(double));
}

S21Matrix::S21Matrix(S21Matrix&& src) : rows_(src.rows_), cols_(src.cols_) {
  //	std::cout << "Constructor: [" << this << "]" << std::endl;
  matrix_len = rows_ * cols_;
  matrix_ = src.matrix_;
  src.matrix_ = nullptr;
  src.rows_ = 0;
  src.cols_ = 0;
}

S21Matrix::~S21Matrix() {
  //	std::cout << "Desctructor: [" << this << "]" << std::endl;
  if (matrix_) delete matrix_;
  matrix_ = nullptr;
}

//=========The constructors below is for testing purposes only!===========

S21Matrix::S21Matrix(Eigen::MatrixXd& src)
    : rows_(src.rows()), cols_(src.cols()) {
  Eigen::MatrixXd m;

  m = src.transpose().eval();
  matrix_len = rows_ * cols_;
  matrix_ = new double[matrix_len]();
  std::memcpy(matrix_, m.data(), matrix_len * sizeof(double));
}

//======================Accessors And Mutators==========================

int S21Matrix::getRows() { return rows_; }

int S21Matrix::getCols() { return cols_; }

void S21Matrix::editRows(int new_rows) {
  int k = 0;
  S21Matrix new_matrix(new_rows, cols_);

  if (rows_ >= new_rows) {
    for (int j = 0; j < matrix_len; j++) {
      if (j / cols_ < new_rows) new_matrix.matrix_[k++] = matrix_[j];
    }
  } else {
    for (int j = 0; j < matrix_len; j++) new_matrix.matrix_[k++] = matrix_[j];
  }
  this->AssignMatrix(new_matrix);
}

void S21Matrix::editCols(int new_cols) {
  int j = 0;
  int k = 0;
  S21Matrix new_matrix(rows_, new_cols);

  if (cols_ >= new_cols) {
    for (; j < matrix_len; j++) {
      if (j % cols_ < new_cols) new_matrix.matrix_[k++] = matrix_[j];
    }
  } else {
    for (; k < new_matrix.matrix_len; k++)
      if (cols_ > k % new_cols) new_matrix.matrix_[k] = matrix_[j++];
  }
  this->AssignMatrix(new_matrix);
}

double S21Matrix::getValueByIndex(int index) { return matrix_[index]; }

double& S21Matrix::setValueByIndex(int index) { return matrix_[index]; }

//======================Utility Methods==========================

int S21Matrix::CustomPower(int init_value, int power) {
  int result = init_value;
  if (power == 0) return 1;
  for (int i = 1; i < power; i++) {
    result *= init_value;
  }
  return result;
}

void S21Matrix::FillMatrix(int limit) {
  std::srand(time(0));
  for (int i = 0; i < matrix_len; i++) matrix_[i] = std::rand() % limit;
}

void S21Matrix::ShowMatrix() {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) std::cout << matrix_[cols_ * i + j] << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

void S21Matrix::AssignMatrix(const S21Matrix& other) {
  if (matrix_) delete matrix_;
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_len = rows_ * cols_;
  matrix_ = new double[matrix_len]();
  std::memcpy(matrix_, other.matrix_, matrix_len * sizeof(double));
}

void S21Matrix::RoundZeroes() {
  for (int i = 0; i < matrix_len; i++) {
    if (matrix_[i] < 0.001) matrix_[i] = 0;
  }
}

void S21Matrix::CheckInverse() {
  this->MulMatrix(this->Inverse());
  this->RoundZeroes();
  this->ShowMatrix();
}

double S21Matrix::Det2x2() {
  return matrix_[0] * matrix_[3] - matrix_[1] * matrix_[2];
}

void S21Matrix::Complements2x2(S21Matrix& minors) {
  minors.matrix_[0] = matrix_[3];
  minors.matrix_[1] = -matrix_[2];
  minors.matrix_[2] = -matrix_[1];
  minors.matrix_[3] = matrix_[0];
}

S21Matrix S21Matrix::CropMatrix(int i) {
  int k = 0;
  int crop_row = i / cols_;
  int crop_col = i % cols_;
  S21Matrix croped_matrix(rows_ - 1, cols_ - 1);

  for (int j = 0; j < matrix_len; j++) {
    if (j / cols_ != crop_row && j % cols_ != crop_col)
      croped_matrix.matrix_[k++] = matrix_[j];
  }
  return croped_matrix;
}

void S21Matrix::Describe() {
  std::cout << "\nInitial matrix" << std::endl;
  this->ShowMatrix();

  std::cout << "\nRows: " << rows_ << " Cols: " << cols_
            << " Len: " << matrix_len << " Det: " << this->Determinant()
            << std::endl;

  std::cout << "\nTranspose of the matrix" << std::endl;
  this->Transpose().ShowMatrix();

  std::cout << "\nAlgebraic addition matrix" << std::endl;
  this->CalcComplements().ShowMatrix();

  std::cout << "\nInverse of matrix" << std::endl;
  this->Inverse().ShowMatrix();

  std::cout << "\nInitial Matrix x Inverse Matrix" << std::endl;
  this->CheckInverse();
}

//============================Methods================================

bool S21Matrix::EqMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument(
        "EqMatrix: Incorrect dimensions for matrix comparation.");
  for (int i = 0; i < matrix_len; i++) {
    if (abs(matrix_[i] - other.matrix_[i]) > EPSILON) return 0;
  }
  return 1;
}

void S21Matrix::SumMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument(
        "SumMatrix: Incorrect dimensions for matrix addition.");
  for (int i = 0; i < matrix_len; i++) matrix_[i] += other.matrix_[i];
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  if (rows_ != other.rows_ || cols_ != other.cols_)
    throw std::invalid_argument(
        "SubMatrix: Incorrect dimensions for matrix subtraction.");
  for (int i = 0; i < matrix_len; i++) matrix_[i] -= other.matrix_[i];
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; i < matrix_len; i++) matrix_[i] *= num;
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  S21Matrix aux_matrix(rows_, other.cols_);

  if (cols_ != other.rows_)
    throw std::invalid_argument(
        "MulMatrix: Incorrect dimensions for matrix multiplication.");
  for (int row_A = 0; row_A < rows_; row_A++) {
    for (int col_B = 0; col_B < other.cols_; col_B++) {
      for (int shared_index = 0; shared_index < cols_; shared_index++)
        aux_matrix.matrix_[aux_matrix.cols_ * row_A + col_B] +=
            matrix_[cols_ * row_A + shared_index] *
            other.matrix_[other.cols_ * shared_index + col_B];
    }
  }
  this->AssignMatrix(aux_matrix);
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix matrix_t(cols_, rows_);

  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++)
      matrix_t.matrix_[rows_ * j + i] = matrix_[cols_ * i + j];
  }
  return matrix_t;
}

double S21Matrix::Determinant() {
  double det = 0;
  if (rows_ == 1 && cols_ == 1) return matrix_[0];
  if (rows_ != cols_)
    throw std::invalid_argument("Determinant: Matrix is not square.");
  if (rows_ >= 3) {
    for (int i = 0; i < cols_; i++)
      det += matrix_[i] * this->CropMatrix(i).Determinant() *
             this->CustomPower(-1, i / cols_ + i % cols_);
  } else
    return this->Det2x2();
  return det;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix minors(rows_, cols_);

  if (rows_ == 1 && cols_ == 1) return *this;
  if (rows_ != cols_)
    throw std::invalid_argument("CalcComplements: Matrix is not square.");
  if (rows_ >= 3) {
    for (int i = 0; i < matrix_len; i++)
      minors.matrix_[i] = this->CropMatrix(i).Determinant() *
                          this->CustomPower(-1, i / cols_ + i % cols_);
  } else
    this->Complements2x2(minors);
  return minors;
}

S21Matrix S21Matrix::Inverse() {
  S21Matrix matrix_i;

  if (rows_ == 1 && cols_ == 1) {
    if (matrix_[0] == 0)
      return *this;
    else {
      matrix_[0] = 1 / matrix_[0];
      return *this;
    }
  }
  double det = this->Determinant();

  if (det == 0)
    throw std::invalid_argument("Inverse: Determinant equals to zero!");
  matrix_i.AssignMatrix(this->CalcComplements().Transpose());
  matrix_i.MulNumber(1 / det);
  return matrix_i;
}

//=========================== Operators =================================

S21Matrix S21Matrix::operator+(const S21Matrix& other) {
  S21Matrix copy_matrix(*this);
  copy_matrix.SumMatrix(other);
  return copy_matrix;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) {
  S21Matrix copy_matrix(*this);
  copy_matrix.SubMatrix(other);
  return copy_matrix;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) {
  S21Matrix copy_matrix(*this);
  copy_matrix.MulMatrix(other);
  return copy_matrix;
}

bool S21Matrix::operator==(const S21Matrix& other) {
  return this->EqMatrix(other);
}

void S21Matrix::operator=(const S21Matrix& other) { this->AssignMatrix(other); }

void S21Matrix::operator+=(const S21Matrix& other) { this->SumMatrix(other); }

void S21Matrix::operator-=(const S21Matrix& other) { this->SubMatrix(other); }

void S21Matrix::operator*=(const S21Matrix& other) { this->MulMatrix(other); }

void S21Matrix::operator*=(double num) { this->MulNumber(num); }

double S21Matrix::operator()(int i, int j) {
  int index = cols_ * i + j;

  if (index >= matrix_len || index < 0)
    throw std::out_of_range("operator(): Index out of range");
  return matrix_[index];
}
