import sympy
from sympy import *
import numpy

sympy.init_printing(use_unicode=True)

x, y, z = sympy.symbols('x y z')

a = sympy.Rational(1, 2)
a * 2
sympy.pi**2

# evalf evaluates the expression to a floating-point number
sympy.pi.evalf()
(sympy.pi + sympy.exp(1)).evalf()
cos(2 * x).evalf(subs={x: 2.4})

# convert a SymPy expression to an expression that can be numerically evaluated
a = numpy.arange(10)
expr = sympy.sin(x)
f = lambdify(x, expr, "numpy")
f(a)