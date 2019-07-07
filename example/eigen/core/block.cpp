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

void block_as_rvalue();
void block_as_lvalue();
void row_column_as_special_block();
void corner_related_block();
void block_for_vector();

int main() {
    block_as_rvalue();
    block_as_lvalue();
    row_column_as_special_block();
    corner_related_block();
    block_for_vector();
}

void block_as_rvalue() {
    Eigen::MatrixXf m(4, 4);
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    cout << "block in the middle" << endl;
    cout << m.block<2, 2>(1, 1) << endl << endl;
    for (int i = 1; i <= 3; ++i) {
        cout << "block of size " << i << "x" << i << endl;
        cout << m.block(0, 0, i, i) << endl << endl;
    }
}

void block_as_lvalue() {
    Array22f m;
    m << 1, 2, 3, 4;
    Array44f a = Array44f::Constant(0.6);
    cout << "Here is the array a:" << endl << a << endl << endl;
    a.block<2, 2>(1, 1) = m;
    cout << "Here is now a with m copied into its central 2x2 block:" << endl << a << endl << endl;
    a.block(0, 0, 2, 3) = a.block(2, 1, 2, 3);
    cout << "Here is now a with bottom-right 2x3 block copied into top-left 2x3 block:" << endl << a << endl << endl;
}

void row_column_as_special_block() {
    Eigen::MatrixXf m(3, 3);
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    cout << "Here is the matrix m:" << endl << m << endl;
    cout << "2nd Row: " << m.row(1) << endl;
    m.col(2) += 3 * m.col(0);
    cout << "After adding 3 times the first column into the third column, the matrix m is:\n";
    cout << m << endl;
}

void corner_related_block() {
    /**
     * matrix.topLeftCorner(p,q);
     * matrix.bottomLeftCorner(p,q);
     * matrix.topRightCorner(p,q);
     * matrix.bottomRightCorner(p,q);
     * matrix.topRows(q);
     * matrix.bottomRows(q);
     * matrix.leftCols(p);
     * matrix.rightCols(q);
     */
    Eigen::Matrix4f m;
    m << 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16;
    cout << "m.leftCols(2) =" << endl << m.leftCols(2) << endl << endl;
    cout << "m.bottomRows<2>() =" << endl << m.bottomRows<2>() << endl << endl;
    m.topLeftCorner(1, 3) = m.bottomRightCorner(3, 1).transpose();
    cout << "After assignment, m = " << endl << m << endl << endl;
}

void block_for_vector() {
    /**
     * block containing the first n elements, vector.head(n);
     * block containing the last n elements, vector.tail(n);
     * block containing n elements, starting at position i, vector.segment(i,n);
     */
    Eigen::ArrayXf v(6);
    v << 1, 2, 3, 4, 5, 6;
    cout << "v.head(3) =" << endl << v.head(3) << endl << endl;
    cout << "v.tail<3>() = " << endl << v.tail<3>() << endl << endl;
    v.segment(1, 4) *= 2;
    cout << "after 'v.segment(1,4) *= 2', v =" << endl << v << endl;
}
