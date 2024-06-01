import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

x, y, z = sympy.symbols('x y z')

# Polynomial/Rational Function Simplification
sympy.expand((x + y) ** 3)
expanded_expr = sympy.expand(x*(x + 2*y))
sympy.expand(sympy.cos(x + y), trig=True)
sympy.factor(expanded_expr)
sympy.simplify((x + x * y) / x)
simplify((x**3 + x**2 - x - 1)/(x**2 + 2*x + 1))
collect(x*y + x - 3 + 2*x**2 - z*x**2 + x**3, x)
cancel((x**2 + 2*x + 1)/(x**2 + x))
cancel((x*y**2 - 2*x*y*z + x*z**2 + y**2 - 2*y*z + z**2)/(x**2 - 1))
apart((4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x)) # partial fraction decomposition on a rational function
