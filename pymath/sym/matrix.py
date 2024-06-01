from sympy import *
import sympy
init_printing(use_unicode=True)

Matrix([[1, -1], [3, 4], [0, 2]])
Matrix([1, 2, 3])  # column matrix
eye(3)
zeros(2, 3)
ones(3, 2)
diag(1, 2, 3)
M = Matrix([[1, 2, 3], [3, 2, 1]])

M.shape
M.row(0) # first row
M.col(-1) # last column

M.col_del(0) # in place delete
M.row_del(1) # in place delete
M = M.row_insert(1, Matrix([[0, 4]])) # insert do not operate in place
M = M.col_insert(0, Matrix([1, -2]))

M = Matrix([[1, 3], [-2, 3]])
N = Matrix([[0, 3], [0, 7]])
M + N
M * N
3*M
M.T # transpose
M**2  # M * M
M**-1  # invert
N**-1  # not invertible
M.det() # determinant
Matrix([[1, 0, 1, 3], [2, 3, 4, 7], [-1, -3, -3, -4]]).rref() # reduced row echelon form
Matrix([[1, 2, 3, 0, 0], [4, 10, 0, 0, 1]]).nullspace() # a list of column vectors that span the nullspace
Matrix([[1, 1, 2], [2 ,1 , 3], [3 , 1, 4]]).columnspace()

# eigen
M = Matrix([[3, -2,  4, -2], [5,  3, -3, -2], [5, -2,  2, -2], [5, -2, -3,  3]])
M.eigenvals()
M.eigenvects()
lamda = symbols('lamda')
p = M.charpoly(lamda)
factor(p.as_expr())

# factor
P, D = M.diagonalize()  # M = PDP^-1
P*D*P**-1 == M
M.LUdecomposition()
M.LDLdecomposition()