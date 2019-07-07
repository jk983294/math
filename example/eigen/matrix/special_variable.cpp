#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * a block is a rectangular part of a matrix or array
 *
 * block of size (p, q), starting at (i, j)
 * matrix.block(i, j, p, q);
 * matrix.block<p, q>(i, j);
 */

void zero();
void constant();
void random_data();
void identity();
void lin_spaced();
void set_special();

int main() {
    zero();
    constant();
    random_data();
    identity();
    lin_spaced();
    set_special();
}

void zero() {
    std::cout << "A fixed-size array:\n";
    Array33f a1 = Array33f::Zero();
    std::cout << a1 << "\n\n";
    std::cout << "A one-dimensional dynamic-size array:\n";
    ArrayXf a2 = ArrayXf::Zero(3);
    std::cout << a2 << "\n\n";
    std::cout << "A two-dimensional dynamic-size array:\n";
    ArrayXXf a3 = ArrayXXf::Zero(3, 4);
    std::cout << a3 << "\n\n";
}

void constant() {
    MatrixXd m = MatrixXd::Constant(2, 2, 3);
    std::cout << "A const matrix:\n" << m << "\n\n";

    ArrayXXf a1 = ArrayXXf::Constant(1, 2, 3);
    std::cout << "A const array:\n" << a1 << "\n\n";
}

/**
 * numbers are uniformly spread through their whole definition range for integer types,
 * and in the [-1:1] range for floating point scalar types
 */
void random_data() {
    MatrixXd m = MatrixXd::Random(2, 3);
    std::cout << "A random matrix:\n" << m << "\n\n";

    ArrayXXf a1 = ArrayXXf::Random(1, 3);
    std::cout << "A random array:\n" << a1 << "\n\n";

    cout << "A random vector:\n" << VectorXi::Random(2) << "\n\n";
}

void identity() { cout << "A identity matrix:\n" << MatrixXd::Identity(4, 3) << "\n\n"; }

/**
 * LinSpaced(size, low, high) is only available for vectors and one-dimensional arrays
 */
void lin_spaced() {
    ArrayXXf table(10, 4);
    table.col(0) = ArrayXf::LinSpaced(10, 0, 90);
    table.col(1) = M_PI / 180 * table.col(0);
    table.col(2) = table.col(1).sin();
    table.col(3) = table.col(1).cos();
    std::cout << "  Degrees   Radians      Sine    Cosine\n";
    std::cout << table << "\n\n";
}

/**
 * three different ways to generate the same matrix
 */
void set_special() {
    const int size = 6;
    MatrixXd mat1(size, size);
    mat1.topLeftCorner(size / 2, size / 2) = MatrixXd::Zero(size / 2, size / 2);
    mat1.topRightCorner(size / 2, size / 2) = MatrixXd::Identity(size / 2, size / 2);
    mat1.bottomLeftCorner(size / 2, size / 2) = MatrixXd::Identity(size / 2, size / 2);
    mat1.bottomRightCorner(size / 2, size / 2) = MatrixXd::Zero(size / 2, size / 2);
    std::cout << mat1 << std::endl << std::endl;

    MatrixXd mat2(size, size);
    mat2.topLeftCorner(size / 2, size / 2).setZero();
    mat2.topRightCorner(size / 2, size / 2).setIdentity();
    mat2.bottomLeftCorner(size / 2, size / 2).setIdentity();
    mat2.bottomRightCorner(size / 2, size / 2).setZero();
    std::cout << mat2 << std::endl << std::endl;

    MatrixXd mat3(size, size);
    mat3 << MatrixXd::Zero(size / 2, size / 2), MatrixXd::Identity(size / 2, size / 2),
        MatrixXd::Identity(size / 2, size / 2), MatrixXd::Zero(size / 2, size / 2);
    std::cout << mat3 << std::endl;
}
