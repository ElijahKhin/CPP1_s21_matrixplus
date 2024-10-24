#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "s21_matrix_oop.h"

using namespace Eigen;

TEST(Constructor, BasicConstr) {
  S21Matrix m;

  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 2);
  EXPECT_EQ(m.matrix_len, 4);
}

TEST(Constructor, ParamConstrException) {
  try {
    S21Matrix m(0, 0);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string("Error in parameterized constructor: ∀rows_,cols_ ∈ "
                          "[1, 46340]\nExit Failure\n"),
              e.what());
  }
  try {
    S21Matrix m(46341, 0);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string("Error in parameterized constructor: ∀rows_,cols_ ∈ "
                          "[1, 46340]\nExit Failure\n"),
              e.what());
  }
}

TEST(Constructor, CopyConstr) {
  MatrixXd e_matrix;
  e_matrix = (MatrixXd::Random(3, 3).array() * 100);

  S21Matrix m(e_matrix);

  S21Matrix m_copy(m);

  EXPECT_TRUE(m == m_copy);
}

TEST(Constructor, MoveConstr) {
  MatrixXd e_matrix;
  e_matrix = (MatrixXd::Random(3, 3).array() * 100);

  S21Matrix m_l(e_matrix);
  S21Matrix m_r(e_matrix);

  S21Matrix m_move(m_l + m_r);

  EXPECT_TRUE(m_move == (m_l + m_r));
}

TEST(Accessors, getRows) {
  S21Matrix m(2, 3);

  EXPECT_EQ(2, m.getRows());
}

TEST(Accessors, getCols) {
  S21Matrix m(2, 3);

  EXPECT_EQ(3, m.getCols());
}

TEST(Mutators, editRows) {
  S21Matrix m;

  m.editRows(3);
  EXPECT_EQ(3, m.getRows());
  m.editRows(1);
  EXPECT_EQ(1, m.getRows());
}

TEST(Mutators, editCols) {
  S21Matrix m;

  m.editCols(3);
  EXPECT_EQ(3, m.getCols());
  m.editCols(1);
  EXPECT_EQ(1, m.getCols());
}

TEST(EqMatrix, Dim1) {
  MatrixXd e_matrix;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix = (MatrixXd::Random(1, 1).array() * 100);
    S21Matrix m1(e_matrix);
    S21Matrix m2(e_matrix);

    EXPECT_TRUE(m1.EqMatrix(m2));

    m2.setValueByIndex(0) = m2.getValueByIndex(0) + 1;
    EXPECT_FALSE(m1.EqMatrix(m2));
  }
}

TEST(EqMatrix, Dim2) {
  MatrixXd e_matrix;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix = (MatrixXd::Random(2, 2).array() * 100);
    S21Matrix m1(e_matrix);
    S21Matrix m2(e_matrix);

    EXPECT_TRUE(m1.EqMatrix(m2));

    m2.setValueByIndex(0) = m2.getValueByIndex(0) + 1;
    EXPECT_FALSE(m1.EqMatrix(m2));
  }
}

TEST(EqMatrix, DimRandom) {
  MatrixXd e_matrix;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    S21Matrix m1(e_matrix);
    S21Matrix m2(e_matrix);

    EXPECT_TRUE(m1.EqMatrix(m2));
    m2.setValueByIndex(0) = m2.getValueByIndex(0) + 1;
    EXPECT_FALSE(m1.EqMatrix(m2));
  }
}

TEST(EqMatrix, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(2, 3);

  try {
    m1.EqMatrix(m2);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(
        std::string("EqMatrix: Incorrect dimensions for matrix comparation."),
        e.what());
  }
}

TEST(SumMatrix, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SumMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SumMatrix, Dim2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_r = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SumMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SumMatrix, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix_l = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SumMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SumMatrix, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(2, 3);

  try {
    m1.SumMatrix(m2);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(
        std::string("SumMatrix: Incorrect dimensions for matrix addition."),
        e.what());
  }
}

TEST(SubMatrix, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SubMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SubMatrix, Dim2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_r = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SubMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SubMatrix, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix_l = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.SubMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(SubMatrix, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(2, 3);

  try {
    m1.SubMatrix(m2);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(
        std::string("SubMatrix: Incorrect dimensions for matrix subtraction."),
        e.what());
  }
}

TEST(MulNumber, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l * i;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_l.MulNumber(i);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulNumber, Dim2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l * i;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_l.MulNumber(i);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulNumber, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix_l = (MatrixXd::Random(rows + 1, rows + 1).array() * 100);
    e_matrix_res = e_matrix_l * i;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_l.MulNumber(i);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulMatrix, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.MulMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulMatrix, Dim2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_r = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.MulMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulMatrix, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_r;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix_l = (MatrixXd::Random(rows + 1, rows + 1).array() * 100);
    e_matrix_r = (MatrixXd::Random(rows + 1, cols + 1).array() * 100);
    e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    m_l.MulMatrix(m_r);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulMatrix, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(3, 3);

  try {
    m1.MulMatrix(m2);
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string(
                  "MulMatrix: Incorrect dimensions for matrix multiplication."),
              e.what());
  }
}

TEST(Transpose, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;
  S21Matrix m_r;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l.transpose().eval();

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_r = m_l.Transpose();

    EXPECT_TRUE(m_r.EqMatrix(m_res));
  }
}

TEST(Transpose, Dim2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;
  S21Matrix m_r;

  std::srand(time(0));
  for (int i = 0; i < 1; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l.transpose().eval();

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_r = m_l.Transpose();

    EXPECT_TRUE(m_r.EqMatrix(m_res));
  }
}

TEST(Transpose, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;
  S21Matrix m_r;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int rows = std::rand() % 100;
    int cols = std::rand() % 100;
    e_matrix_l = (MatrixXd::Random(rows + 1, rows + 1).array() * 100);
    e_matrix_res = e_matrix_l.transpose().eval();

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    m_r = m_l.Transpose();

    EXPECT_TRUE(m_r.EqMatrix(m_res));
  }
}

TEST(Determinant, Dim1) {
  MatrixXd e_matrix_l;
  double e_matrix_res;
  double m_r;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 100);
    e_matrix_res = e_matrix_l.determinant();

    S21Matrix m_l(e_matrix_l);

    m_r = m_l.Determinant();

    EXPECT_TRUE(abs(m_r - e_matrix_res) < EPSILON);
  }
}

TEST(Determinant, Dim2) {
  MatrixXd e_matrix_l;
  double e_matrix_res;
  double m_r;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l.determinant();

    S21Matrix m_l(e_matrix_l);

    m_r = m_l.Determinant();

    EXPECT_TRUE(abs(m_r - e_matrix_res) < EPSILON);
  }
}

TEST(Determinant, DimRandom) {
  MatrixXd e_matrix_l;
  double e_matrix_res;
  double m_r;

  std::srand(time(0));
  for (int i = 0; i < 1000; i++) {
    int dim = std::rand() % 8;
    e_matrix_l = (MatrixXd::Random(dim + 1, dim + 1).array() * 10);
    e_matrix_res = e_matrix_l.determinant();

    S21Matrix m_l(e_matrix_l);

    m_r = m_l.Determinant();

    EXPECT_TRUE(abs(m_r - e_matrix_res) < EPSILON);
  }
}

TEST(Determinant, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(2, 3);

  try {
    m1.Determinant();
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string("Determinant: Matrix is not square."), e.what());
  }
}

TEST(Inverse, ZeroCheck) {
  MatrixXd e_matrix_l(1, 1);
  MatrixXd e_matrix_res;
  S21Matrix m_r;

  e_matrix_l << 0;
  e_matrix_res = e_matrix_l.inverse();
  S21Matrix m_l(e_matrix_l);
  m_r = m_l.Inverse();
  EXPECT_TRUE(m_r.EqMatrix(m_l));
}

TEST(Inverse, Dim1) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;
  S21Matrix m_r;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    e_matrix_res = e_matrix_l.inverse();

    S21Matrix m_l(e_matrix_l);

    m_r = m_l.Inverse();

    EXPECT_TRUE(m_r.EqMatrix(m_l));
  }
}

TEST(Inverse, Dim2x2) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    e_matrix_l = (MatrixXd::Random(2, 2).array() * 100);
    e_matrix_res = e_matrix_l.inverse();
    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);
    S21Matrix m_r(m_l.Inverse());

    EXPECT_TRUE(m_r.EqMatrix(m_res));
  }
}

TEST(Inverse, DimRandom) {
  MatrixXd e_matrix_l;
  MatrixXd e_matrix_res;

  std::srand(time(0));
  for (int i = 0; i < 10; i++) {
    int dim = std::rand() % 8;

    e_matrix_l = (MatrixXd::Random(dim + 1, dim + 1).array() * 100);
    e_matrix_res = e_matrix_l.inverse();
    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);
    S21Matrix m_r(m_l.Inverse());

    EXPECT_TRUE(m_r.EqMatrix(m_res));
  }
}

TEST(Inverse, Exception) {
  S21Matrix m1(3, 3);

  try {
    m1.Inverse();
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string("Inverse: Determinant equals to zero!"), e.what());
  }
}

TEST(CalcComplements, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(m_l.CalcComplements());

    EXPECT_EQ(m_r.getValueByIndex(0), e_matrix_l(0, 0));
  }
}

TEST(CalcComplements, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(m_l.CalcComplements());

    EXPECT_EQ(m_r.getValueByIndex(0), e_matrix_l(1, 1));
    EXPECT_EQ(m_r.getValueByIndex(1), -e_matrix_l(1, 0));
    EXPECT_EQ(m_r.getValueByIndex(2), -e_matrix_l(0, 1));
    EXPECT_EQ(m_r.getValueByIndex(3), e_matrix_l(0, 0));
  }
}

TEST(CalcComplements, Dim3x3) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(3, 3).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(m_l.CalcComplements());

    EXPECT_EQ(m_r.getValueByIndex(0), e_matrix_l({1, 2}, {1, 2}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(1),
              e_matrix_l({1, 2}, {0, 2}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(2), e_matrix_l({1, 2}, {0, 1}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(3),
              e_matrix_l({0, 2}, {1, 2}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(4), e_matrix_l({0, 2}, {0, 2}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(5),
              e_matrix_l({0, 2}, {0, 1}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(6), e_matrix_l({0, 1}, {1, 2}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(7),
              e_matrix_l({0, 1}, {0, 2}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(8), e_matrix_l({0, 1}, {0, 1}).determinant());
  }
}

TEST(CalcComplements, Dim4x4) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(4, 4).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(m_l.CalcComplements());

    EXPECT_EQ(m_r.getValueByIndex(0),
              e_matrix_l({1, 2, 3}, {1, 2, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(1),
              e_matrix_l({1, 2, 3}, {0, 2, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(2),
              e_matrix_l({1, 2, 3}, {0, 1, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(3),
              e_matrix_l({1, 2, 3}, {0, 1, 2}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(4),
              e_matrix_l({0, 2, 3}, {1, 2, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(5),
              e_matrix_l({0, 2, 3}, {0, 2, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(6),
              e_matrix_l({0, 2, 3}, {0, 1, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(7),
              e_matrix_l({0, 2, 3}, {0, 1, 2}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(8),
              e_matrix_l({0, 1, 3}, {1, 2, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(9),
              e_matrix_l({0, 1, 3}, {0, 2, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(10),
              e_matrix_l({0, 1, 3}, {0, 1, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(11),
              e_matrix_l({0, 1, 3}, {0, 1, 2}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(12),
              e_matrix_l({0, 1, 2}, {1, 2, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(13),
              e_matrix_l({0, 1, 2}, {0, 2, 3}).determinant());
    EXPECT_EQ(m_r.getValueByIndex(14),
              e_matrix_l({0, 1, 2}, {0, 1, 3}).determinant() * -1);
    EXPECT_EQ(m_r.getValueByIndex(15),
              e_matrix_l({0, 1, 2}, {0, 1, 2}).determinant());
  }
}

TEST(CalcComplements, Exception) {
  S21Matrix m1(3, 2);
  S21Matrix m2(2, 3);

  try {
    m1.CalcComplements();
  } catch (std::invalid_argument& e) {
    EXPECT_EQ(std::string("CalcComplements: Matrix is not square."), e.what());
  }
}

TEST(SumOperator, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l + m_r).EqMatrix(m_res));
  }
}

TEST(SumOperator, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l + m_r).EqMatrix(m_res));
  }
}

TEST(SumOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l + e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l + m_r).EqMatrix(m_res));
  }
}

TEST(SubOperator, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l - m_r).EqMatrix(m_res));
  }
}

TEST(SubOperator, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l - m_r).EqMatrix(m_res));
  }
}

TEST(SubOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l - e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l - m_r).EqMatrix(m_res));
  }
}

TEST(MulMatrixOperator, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * m_r).EqMatrix(m_res));
  }
}

TEST(MulMatrixOperator, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * m_r).EqMatrix(m_res));
  }
}

TEST(MulMatrixOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim2 + 1, dim1 + 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * e_matrix_r;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * m_r).EqMatrix(m_res));
  }
}

TEST(MulNumberRightOperator, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * 5;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * 5).EqMatrix(m_res));
  }
}

TEST(MulNumberRightOperator, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * 5;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * 5).EqMatrix(m_res));
  }
}

TEST(MulNumberRightOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_res = e_matrix_l * 5;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((m_l * 5).EqMatrix(m_res));
  }
}

TEST(MulNumberLeftOperator, Dim1) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(1, 1).array() * 10);
    MatrixXd e_matrix_res = 5 * e_matrix_l;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((5 * m_l).EqMatrix(m_res));
  }
}

TEST(MulNumberLeftOperator, Dim2x2) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    MatrixXd e_matrix_l = (MatrixXd::Random(2, 2).array() * 10);
    MatrixXd e_matrix_res = 5 * e_matrix_l;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((5 * m_l).EqMatrix(m_res));
  }
}

TEST(MulNumberLeftOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 100; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_res = 5 * e_matrix_l;

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_res(e_matrix_res);

    EXPECT_TRUE((5 * m_l).EqMatrix(m_res));
  }
}

TEST(EqOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);

    EXPECT_EQ(m_l == m_r, e_matrix_l == e_matrix_r);
  }
}

TEST(AssignmentOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);

    m_l = m_r;
    e_matrix_l = e_matrix_r;

    S21Matrix m_res(e_matrix_l);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(AddAssignmentOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);

    m_l += m_r;
    e_matrix_l += e_matrix_r;

    S21Matrix m_res(e_matrix_l);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(DiffAssignmentOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);

    m_l -= m_r;
    e_matrix_l -= e_matrix_r;

    S21Matrix m_res(e_matrix_l);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulMatrixAssignmentOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);
    MatrixXd e_matrix_r = (MatrixXd::Random(dim2 + 1, dim1 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    S21Matrix m_r(e_matrix_r);

    m_l *= m_r;
    e_matrix_l *= e_matrix_r;

    S21Matrix m_res(e_matrix_l);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(MulNumberAssignmentOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);

    m_l *= 10;
    e_matrix_l *= 10;

    S21Matrix m_res(e_matrix_l);

    EXPECT_TRUE(m_l.EqMatrix(m_res));
  }
}

TEST(IndexationOperator, DimRandom) {
  std::srand(time(0));
  for (int i = 0; i < 3000; i++) {
    int dim1 = std::rand() % 100;
    int dim2 = std::rand() % 100;
    MatrixXd e_matrix_l = (MatrixXd::Random(dim1 + 1, dim2 + 1).array() * 10);

    S21Matrix m_l(e_matrix_l);
    for (int i = 0; i <= dim1; i++) {
      for (int j = 0; j <= dim2; j++) EXPECT_EQ(m_l(i, j), e_matrix_l(i, j));
    }
  }
}

TEST(IndexationOperator, Exception) {
  S21Matrix m_l;

  try {
    m_l(-1, -1);
  } catch (std::out_of_range& e) {
    EXPECT_EQ(std::string("operator(): Index out of range"), e.what());
  }
  try {
    m_l(MAX_DIM + 1, MAX_DIM + 1);
  } catch (std::out_of_range& e) {
    EXPECT_EQ(std::string("operator(): Index out of range"), e.what());
  }
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
