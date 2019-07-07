#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * a matrix is stored in row-major order if it is stored row by row
 * a matrix is stored in column-major order if it is stored column by column
 *
 * Eigen defaults to storing the entry in column-major
 */

int main() {
    Matrix<int, 3, 4, ColMajor> colMajorMatrix;
    colMajorMatrix << 8, 2, 2, 9, 9, 1, 4, 4, 3, 5, 4, 5;
    cout << "The matrix A:" << endl;
    cout << colMajorMatrix << endl << endl;
    cout << "In memory (column-major):" << endl;
    for (int i = 0; i < colMajorMatrix.size(); i++) cout << *(colMajorMatrix.data() + i) << "  ";
    cout << endl << endl;  // 8  9  3  2  1  5  2  4  4  9  4  5

    Matrix<int, 3, 4, RowMajor> rowMajorMatrix = colMajorMatrix;
    cout << "In memory (row-major):" << endl;
    for (int i = 0; i < rowMajorMatrix.size(); i++) cout << *(rowMajorMatrix.data() + i) << "  ";
    cout << endl;  // 8  2  2  9  9  1  4  4  3  5  4  5
}
