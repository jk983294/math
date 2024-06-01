import sympy
from sympy import *

sympy.init_printing(use_unicode=True)
x, y, z = symbols('x y z')

# solving algebraic equations
sympy.solve(x**2 - 2, x)
solveset(x**2 - x, x)
solveset(x - x, x, domain=S.Reals)
solveset(sin(x) - 1, x, domain=S.Reals)
solveset(exp(x), x)     # No solution exists
solveset(cos(x) - x, x)  # Not able to find solution

# the linear system of equations
linsolve([x + y + z - 1, x + y + 2*z - 3 ], (x, y, z)) # List of Equations Form
linsolve(Matrix(([1, 1, 1, 1], [1, 1, 2, 3])), (x, y, z)) # Augmented Matrix Form
M = Matrix(((1, 1, 1, 1), (1, 1, 2, 3)))
system = A, b = M[:, :-1], M[:, -1]
linsolve(system, x, y, z)  # A*x = b Form

# non linear system of equations
a, b, c, d = symbols('a, b, c, d', real=True)
nonlinsolve([a**2 + a, a - b], [a, b])  # only real solution is present
nonlinsolve([x**2 + 1, y**2 + 1], [x, y])  # only complex solution is present
nonlinsolve([x**2 - 2*y**2 -2, x*y - 2], [x, y])  # both real and complex solution are present
nonlinsolve([x*y, x*y - x], [x, y])  # infinitely many solutions

# solve differential equations
f, g = symbols('f g', cls=Function)
diffeq = Eq(f(x).diff(x, x) - 2*f(x).diff(x) + f(x), sin(x))
dsolve(diffeq, f(x))