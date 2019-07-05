#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * the Array class provides an easy way to perform coefficient-wise operations,
 * which might not have a linear algebraic meaning, such as adding a constant to every coefficient in the array
 * or multiplying two arrays coefficient-wise
 */

void addition();
void multiplication();
void coefficient_wise_operations();
void array_matrix_convert1();
void array_matrix_convert2();

int main() {
    addition();
    multiplication();
    coefficient_wise_operations();
    array_matrix_convert1();
    array_matrix_convert2();
}

void addition() {
    ArrayXXf a(3, 3);
    ArrayXXf b(3, 3);
    a << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    b << 1, 2, 3, 1, 2, 3, 1, 2, 3;
    cout << "a + b = " << endl << a + b << endl << endl;
    // subtracting a scalar from an array
    cout << "a - 2 = " << endl << a - 2 << endl << endl;
}

void multiplication() {
    ArrayXXf a(2, 2);
    ArrayXXf b(2, 2);
    a << 1, 2, 3, 4;
    b << 5, 6, 7, 8;
    // arrays interpret multiplication as coefficient-wise product
    cout << "a * b = " << endl << a * b << endl << endl;
}

void coefficient_wise_operations() {
    ArrayXf a = ArrayXf::Random(5);
    a *= 2;
    cout << "a =" << endl << a << endl;
    cout << "a.abs() =" << endl << a.abs() << endl;
    cout << "a.abs().sqrt() =" << endl << a.abs().sqrt() << endl;
    cout << "a.min(a.abs().sqrt()) =" << endl << a.min(a.abs().sqrt()) << endl << endl;
}

void array_matrix_convert1() {
    MatrixXf m(2, 2);
    MatrixXf n(2, 2);
    MatrixXf result(2, 2);
    m << 1, 2, 3, 4;
    n << 5, 6, 7, 8;
    result = m * n;
    cout << "-- Matrix m*n: --" << endl << result << endl << endl;
    result = m.array() * n.array();
    cout << "-- Array m*n: --" << endl << result << endl << endl;
    result = m.cwiseProduct(n);
    cout << "-- With cwiseProduct: --" << endl << result << endl << endl;
    result = m.array() + 4;
    cout << "-- Array m + 4: --" << endl << result << endl << endl;
}

void array_matrix_convert2() {
    MatrixXf m(2, 2);
    MatrixXf n(2, 2);
    MatrixXf result(2, 2);
    m << 1, 2, 3, 4;
    n << 5, 6, 7, 8;

    result = (m.array() + 4).matrix() * m;
    cout << "-- Combination 1: --" << endl << result << endl << endl;
    result = (m.array() * n.array()).matrix() * m;
    cout << "-- Combination 2: --" << endl << result << endl << endl;
}
