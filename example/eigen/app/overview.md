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