### reference

http://eigen.tuxfamily.org/dox/group__TopicLinearAlgebraDecompositions.html

### solver

| Decomposition        | Requirements   |  Speed(small to median)  |  Speed(large)  |  Accuracy  |
| --------   | -----:  | :----:  | :----:  | :----:  |
| partialPivLu()     | Invertible |   2     |   2     |   1     |
| fullPivLu()        |   None   |   -1  |   -2     |   3     |
| householderQr()        |    None    |  2  |   2     |   1     |
| colPivHouseholderQr()        |    None    |  1  |   -1     |   3     |
| fullPivHouseholderQr()        |    None    |  -1  |   -2     |   3     |
| completeOrthogonalDecomposition()        |    None    |  1  |   -1     |   3     |
| llt()        |    Positive definite    |  3  |   3     |   1     |
| ldlt()        |    Positive or negative semidefinite    |  3  |   1     |   2     |
| bdcSvd()        |    None    |  -1  |   -1     |   3     |
| jacobiSvd()        |    None    |  -1  |   -3     |   3     |

### Terminology
self-adjoint
* For a real matrix, self adjoint is a synonym for symmetric. For a complex matrix, self adjoint is a synonym for hermitian. 
More generally, a matrix A is self adjoint if and only if it is equal to its adjoint A*. The adjoint is also called the conjugate transpose. 

positive/negative definite
* A selfadjoint matrix A is positive definite if v* Av > 0 for any non zero vector v. 
In the same vein, it is negative definite if v* Av < 0  for any non zero vector v

positive/negative semidefinite
* A self adjoint matrix A is positive semi-definite if v* Av >= 0  for any non zero vector v. 
In the same vein, it is negative semi-definite if v^* Av <= 0 for any non zero vector v