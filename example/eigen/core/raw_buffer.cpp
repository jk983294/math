#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * how to work with raw C/C++ arrays
 * particularly useful when importing vectors and matrices from other libraries into Eigen
 */

void how_to_map();
void how_to_use();
void change_underlying_mapped_array();

int main() {
    how_to_map();
    how_to_use();
    change_underlying_mapped_array();
}

void how_to_map() {
    int array[8];
    for (int i = 0; i < 8; ++i) array[i] = i;
    cout << "Column-major:\n" << Map<Matrix<int, 2, 4> >(array) << endl;
    cout << "Row-major:\n" << Map<Matrix<int, 2, 4, RowMajor> >(array) << endl;
    cout << "Row-major using stride:\n" << Map<Matrix<int, 2, 4>, Unaligned, Stride<1, 4> >(array) << "\n\n";
}

void how_to_use() {
    typedef Matrix<float, 1, Dynamic> MatrixType;
    typedef Map<MatrixType> MapType;
    typedef Map<const MatrixType> MapTypeConst;  // a read-only map
    const int n_dims = 5;

    MatrixType m1(n_dims), m2(n_dims);
    m1.setRandom();
    m2.setRandom();
    float *p = &m2(0);                      // get the address storing the data for m2
    MapType m2map(p, m2.size());            // m2map shares data with m2
    MapTypeConst m2mapconst(p, m2.size());  // a read-only accessor for m2
    cout << "m1: " << m1 << endl;
    cout << "m2: " << m2 << endl;
    cout << "Squared euclidean distance: " << (m1 - m2).squaredNorm() << endl;
    cout << "Squared euclidean distance, using map: " << (m1 - m2map).squaredNorm() << endl;
    m2map(3) = 7;  // this will change m2, since they share the same array
    cout << "Updated m2: " << m2 << endl;
    cout << "m2 coefficient 2, constant accessor: " << m2mapconst(2) << "\n\n";
}

/**
 * using the C++ "placement new" syntax
 */
void change_underlying_mapped_array() {
    int data1[] = {1, 2, 3, 4};
    int data2[] = {5, 6, 7, 8, 9};
    Map<RowVectorXi> v(data1, 4);
    cout << "The mapped vector v is: " << v << "\n";
    new (&v) Map<RowVectorXi>(data2, 5);
    cout << "Now v is: " << v << "\n";
}
