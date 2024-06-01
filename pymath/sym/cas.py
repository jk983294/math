import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

# computer algebra systems
x, y, z = sympy.symbols('x y z')
expr = x + 2*y
expr = expr - 1 - x
expr = sympy.sympify("x**2 + 3*x - 1/2") # build from string

# Substitution variable
expr = x + 2*y
expr.subs(x, 2) # use x = 2 to evaluate expr
expr.subs(x, y + 3)
expr = x**3 + 4*x*y - z
expr.subs([(x, 2), (y, 4), (z, 0)])

# print in latex format
sympy.latex(Integral(cos(x)**2, (x, 0, pi)))
pprint(Integral(sqrt(1/x), x), use_unicode=True)