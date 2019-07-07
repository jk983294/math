#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * reduction operations to reduce a given matrix or vector to a single value
 */

void matrix_reduction();
void norm();
void boolean_reduction();
void visitor();  // find max/min position

/**
 * partial reduction can operate column or row wise on a Matrix or Array
 * column-wise operations return a row vector, while row-wise operations return a column vector
 */
void partial_reduction();
void broadcasting();  // apply vector to row/column of matrix or array
void nearest_neighbour();

int main() {
    matrix_reduction();
    norm();
    boolean_reduction();
    visitor();
    partial_reduction();
    broadcasting();
    nearest_neighbour();
}

void matrix_reduction() {
    Eigen::Matrix2d mat;
    mat << 1, 2, 3, 4;
    cout << "Here is mat.sum():       " << mat.sum() << endl;
    cout << "Here is mat.prod():      " << mat.prod() << endl;
    cout << "Here is mat.mean():      " << mat.mean() << endl;
    cout << "Here is mat.minCoeff():  " << mat.minCoeff() << endl;
    cout << "Here is mat.maxCoeff():  " << mat.maxCoeff() << endl;
    cout << "Here is mat.trace():     " << mat.trace() << endl;
    cout << "another way to calculate mat.trace():     " << mat.diagonal().sum() << endl;
}

/**
 * squaredNorm = sigma (x[i] * x[i])
 * norm = sqrt(squaredNorm)
 * lpNorm<1> = sigma( abs(x[i]) )
 * lpNorm<Infinity> = max ( abs(x[i]) )
 */
void norm() {
    VectorXf v(2);
    v << -1, 2;
    cout << "v.squaredNorm() = " << v.squaredNorm() << endl;
    cout << "v.norm() = " << v.norm() << endl;
    cout << "v.lpNorm<1>() = " << v.lpNorm<1>() << endl;
    cout << "v.lpNorm<Infinity>() = " << v.lpNorm<Infinity>() << endl << endl;

    MatrixXf m(2, 2);
    m << 1, -2, -3, 4;
    cout << "m.squaredNorm() = " << m.squaredNorm() << endl;
    cout << "m.norm() = " << m.norm() << endl;
    cout << "m.lpNorm<1>() = " << m.lpNorm<1>() << endl;
    cout << "m.lpNorm<Infinity>() = " << m.lpNorm<Infinity>() << endl << endl;
}

/**
 * all() returns true if all of the coefficients in a given Matrix or Array evaluate to true
 * any() returns true if at least one of the coefficients in a given Matrix or Array evaluates to true
 * count() returns the number of coefficients in a given Matrix or Array that evaluate to true
 */
void boolean_reduction() {
    ArrayXXf a(2, 2);
    a << 1, 2, 3, 4;
    cout << "(a > 0).all()   = " << (a > 0).all() << endl;            // true
    cout << "(a > 0).any()   = " << (a > 0).any() << endl;            // true
    cout << "(a > 0).count() = " << (a > 0).count() << endl << endl;  // 4
    cout << "(a > 2).all()   = " << (a > 2).all() << endl;            // false
    cout << "(a > 2).any()   = " << (a > 2).any() << endl;            // true
    cout << "(a > 2).count() = " << (a > 2).count() << endl << endl;  // 2
}

void visitor() {
    Eigen::MatrixXf m(2, 2);
    m << 1, 2, 3, 4;
    // get location of maximum
    MatrixXf::Index maxRow, maxCol;
    float max = m.maxCoeff(&maxRow, &maxCol);
    // get location of minimum
    MatrixXf::Index minRow, minCol;
    float min = m.minCoeff(&minRow, &minCol);
    cout << "Max: " << max << ", at: " << maxRow << "," << maxCol << endl;
    cout << "Min: " << min << ", at: " << minRow << "," << minCol << endl << endl;
}

void partial_reduction() {
    Eigen::MatrixXf mat(2, 4);
    mat << 1, 2, 6, 9, 3, 1, 7, 2;

    std::cout << "Column's maximum: " << std::endl << mat.colwise().maxCoeff() << std::endl;
    std::cout << "Row's maximum: " << std::endl << mat.rowwise().maxCoeff() << std::endl;

    MatrixXf::Index maxIndex;
    float maxColumnSum = mat.colwise().sum().maxCoeff(&maxIndex);

    std::cout << "Maximum sum at column " << maxIndex << std::endl;
    std::cout << "The corresponding column is: " << std::endl;
    std::cout << mat.col(maxIndex) << std::endl;
    std::cout << "And its sum is is: " << maxColumnSum << std::endl << endl;
}

void broadcasting() {
    Eigen::MatrixXf mat(2, 4);
    Eigen::VectorXf v(2);
    mat << 1, 2, 6, 9, 3, 1, 7, 2;
    v << -1, 1;

    // add v to each column of m
    mat.colwise() += v;

    std::cout << "Broadcasting result: " << std::endl;
    std::cout << mat << std::endl;

    Eigen::VectorXf v1(4);
    v1 << 0, 1, 2, 3;

    // add v to each row of m
    mat.rowwise() += v1.transpose();

    std::cout << "Broadcasting result: " << std::endl;
    std::cout << mat << std::endl << std::endl;
}

void nearest_neighbour() {
    Eigen::MatrixXf m(2, 4);
    Eigen::VectorXf v(2);

    m << 1, 23, 6, 9, 3, 11, 7, 2;
    v << 2, 3;
    MatrixXf::Index index;
    (m.colwise() - v).colwise().squaredNorm().minCoeff(&index);
    cout << "Nearest neighbour is column " << index << ":" << endl;
    cout << m.col(index) << endl;
}
