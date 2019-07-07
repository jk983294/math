#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * multiply two matrices, Eigen assumes that aliasing occurs. If you know there is no aliasing, you can use noalias()
 * In all other situations, Eigen assumes that there is no aliasing issue and thus gives the wrong result
 * if aliasing does in fact occur. To prevent this, you have to use eval() or one of the xxxInPlace() functions.
 *
 * use Inplace method to avoid alias problem
 * MatrixBase::adjoint()	MatrixBase::adjointInPlace()
 * DenseBase::reverse()	    DenseBase::reverseInPlace()
 * LDLT::solve()	        LDLT::solveInPlace()
 * LLT::solve()	            LLT::solveInPlace()
 * TriangularView::solve()	TriangularView::solveInPlace()
 * DenseBase::transpose()	DenseBase::transposeInPlace()
 */

void alias_bug();
void resolve_alias_issue();
void no_need_eval_for_element_wise_operation();
void noalias_example();

int main() {
    alias_bug();
    resolve_alias_issue();
    no_need_eval_for_element_wise_operation();
    noalias_example();
}

void alias_bug() {
    MatrixXi mat(3, 3);
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    cout << "Here is the matrix mat:\n" << mat << endl;

    /**
     * This assignment shows the aliasing problem, mat[2,2] should = 5, but it = 1
     * The problem is that Eigen uses lazy evaluation for mat.topLeftCorner(2,2)
     */
    mat.bottomRightCorner(2, 2) = mat.topLeftCorner(2, 2);
    cout << "After the assignment, mat = \n" << mat << endl;
}

/**
 * Eigen has to evaluate the right-hand side fully into a temporary matrix/array
 * and then assign it to the left-hand side.
 * The function eval() does precisely that.
 */
void resolve_alias_issue() {
    MatrixXi mat(3, 3);
    mat << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    cout << "Here is the matrix mat:\n" << mat << endl;
    mat.bottomRightCorner(2, 2) = mat.topLeftCorner(2, 2).eval();
    cout << "After the assignment, mat = \n" << mat << endl;
}

void no_need_eval_for_element_wise_operation() {
    MatrixXf mat(2, 2);
    mat << 1, 2, 4, 7;
    cout << "Here is the matrix mat:\n" << mat << endl << endl;
    mat = 2 * mat;
    cout << "After 'mat = 2 * mat', mat = \n" << mat << endl << endl;
    mat = mat - MatrixXf::Identity(2, 2);
    cout << "After the subtraction, it becomes\n" << mat << endl << endl;
    ArrayXXf arr = mat;
    arr = arr.square();
    cout << "After squaring, it becomes\n" << arr << endl << endl;

    // Combining all operations in one statement:
    mat << 1, 2, 4, 7;
    mat = (2 * mat - MatrixXf::Identity(2, 2)).array().square();
    cout << "Doing everything at once yields\n" << mat << endl << endl;
}

void noalias_example() {
    MatrixXf matA(2, 2);
    matA << 2, 0, 0, 2;
    matA = matA * matA;
    cout << "already handled alias\n" << matA << endl << endl;

    MatrixXf matB(2, 2);
    matA << 2, 0, 0, 2;

    matB = matA * matA;
    cout << "Simple but not efficient way:\n" << matB << endl << endl;

    matB.noalias() = matA * matA;
    cout << "complicated but efficient way:\n" << matB << endl << endl;
}
