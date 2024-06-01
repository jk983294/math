import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

y = Function('y')
x, t, z, nu = sympy.symbols('x t z nu')
sympy.dsolve(Eq(y(t).diff(t, t) - y(t), exp(t)), y(t))

f, g = symbols('f g', cls=Function)
de = f(x).diff(x, x) + f(x)
dsolve(de, f(x))