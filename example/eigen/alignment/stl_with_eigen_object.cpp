#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <iostream>
#include <map>

using namespace Eigen;
using namespace std;

/**
 * A 16-byte-aligned allocator must be used. Eigen does provide one ready for use: aligned_allocator
 * If you want to use the std::vector container, you need to #include <Eigen/StdVector>
 *
 * issues arise only with fixed-size vectorizable Eigen types and structures having such Eigen objects as member
 * For other Eigen types, such as Vector3f or MatrixXd, no special care is needed when using STL containers
 */

using Int2Vector4fMap =
    std::map<int, Eigen::Vector4f, std::less<int>, Eigen::aligned_allocator<std::pair<const int, Eigen::Vector4f>>>;

using Vector4V4f = std::vector<Eigen::Vector4f, Eigen::aligned_allocator<Eigen::Vector4f>>;

int main() {
    Int2Vector4fMap m;
    auto& v4 = m[1];
    v4 << 1, 2, 3, 4;

    cout << m[1] << endl << endl;

    Vector4V4f v;
    v.push_back(m[1]);
    cout << v[0] << endl << endl;
}
