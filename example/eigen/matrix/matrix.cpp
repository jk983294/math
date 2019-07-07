#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * typedef Matrix<float, 4, 4> Matrix4f;                // 4 * 4
 * typedef Matrix<float, 3, 1> Vector3f;                // 3 * 1
 * typedef Matrix<int, 1, 2> RowVector2i;               // 1 * 2
 * typedef Matrix<double, Dynamic, Dynamic> MatrixXd;   // m * n
 * typedef Matrix<int, Dynamic, 1> VectorXi;            // m * 1
 */

void init();
void addition();
void scalar_multiplication();
void transpose();
void multiplication();
void product();

int main() {
    init();
    addition();
    scalar_multiplication();
    transpose();
    multiplication();
    product();
}

void init() {
    MatrixXf matA(2, 2);
    matA << 1, 2, 3, 4;
    MatrixXf matB(4, 4);
    matB << matA, matA / 10, matA / 10, matA;
    std::cout << "comma init from 4 sub matrix:\n" << matB << std::endl;

    Matrix3f m1;
    m1.row(0) << 1, 2, 3;
    m1.block(1, 0, 2, 2) << 4, 5, 7, 8;
    m1.col(2).tail(2) << 6, 9;
    std::cout << "init from 3 pieces:\n" << m1 << std::endl;
}

void addition() {
    Matrix2d a;
    a << 1, 2, 3, 4;
    MatrixXd b(2, 2);
    b << 2, 3, 1, 4;

    std::cout << "addition and subtraction" << std::endl;
    std::cout << "a + b =\n" << a + b << std::endl;
    std::cout << "a - b =\n" << a - b << std::endl;
    a += b;
    std::cout << "after a += b; now a =\n" << a << std::endl;
    Vector3d v(1, 2, 3);
    Vector3d w(1, 0, 0);
    std::cout << "-v + w - v =\n" << -v + w - v << std::endl;
}

void scalar_multiplication() {
    Matrix2d a;
    a << 1, 2, 3, 4;
    Vector3d v(1, 2, 3);
    std::cout << "\nscalar multiplication and division" << std::endl;
    std::cout << "a * 2.5 =\n" << a * 2.5 << std::endl;
    std::cout << "0.1 * v =\n" << 0.1 * v << std::endl;
    v *= 2;
    std::cout << "after v *= 2; now v =\n" << v << std::endl;
}

void transpose() {
    Matrix2d a;
    a << 1, 2, 3, 4;
    std::cout << "\ntranspose" << std::endl;
    cout << "a^T non in place transpose:\n" << a.transpose() << endl;
    a.transposeInPlace();
    cout << "a^T in place transpose:\n" << a << endl;
}

void multiplication() {
    std::cout << "\nmatrix-matrix and matrix-vector multiplication" << std::endl;
    Matrix2d mat;
    mat << 1, 2, 3, 4;
    Vector2d u(-1, 1), v(2, 0);
    std::cout << "Here is mat*mat:\n" << mat * mat << std::endl;
    std::cout << "Here is mat*u:\n" << mat * u << std::endl;
    std::cout << "Here is u^T*mat:\n" << u.transpose() * mat << std::endl;
    std::cout << "Here is u^T*v:\n" << u.transpose() * v << std::endl;
    std::cout << "Here is u*v^T:\n" << u * v.transpose() << std::endl;
    mat = mat * mat;
    std::cout << "Let's multiply mat by itself, Now mat is mat:\n" << mat << std::endl;
}

void product() {
    std::cout << "\ndot product and cross product" << std::endl;
    Vector3d v(1, 2, 3);
    Vector3d w(0, 1, 2);
    cout << "Dot product: " << v.dot(w) << endl;
    cout << "Cross product:\n" << v.cross(w) << endl;
}
