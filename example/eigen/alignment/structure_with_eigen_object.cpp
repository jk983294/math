#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * If you define a structure having members of fixed-size vectorizable Eigen types,
 * you must overload its "operator new" so that it generates 16-bytes-aligned pointers
 */

class Foo {
public:
    double x;
    Eigen::Vector2d v;  // no need to be first member

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

/**
 * exmaple to show how to conditionnally align, depending on template parameters
 */
template <int n>
class Bar {
public:
    typedef Eigen::Matrix<float, n, 1> Vector;
    enum { NeedsToAlign = (sizeof(Vector) % 16) == 0 };
    Vector v;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW_IF(NeedsToAlign)
};

struct NoAlignObject {
    Eigen::Matrix<double, 2, 1, Eigen::DontAlign> v;  // disable vectorization
};

int main() {
    Foo *foo = new Foo;
    foo->v << 3, 4;
    cout << foo->v << endl << endl;

    Bar<4> *bar4 = new Bar<4>;  // bar4 is guaranteed to be 128bit-aligned
    Bar<3> *bar3 = new Bar<3>;  // bar3 has only the system default alignment guarantee
    bar4->v << 1, 2, 3, 4;
    cout << bar4->v << endl << endl;
    bar3->v << 1, 2, 3;
    cout << bar3->v << endl << endl;

    NoAlignObject *var = new NoAlignObject;
    var->v << 3, 4;
    cout << var->v << endl << endl;
}
