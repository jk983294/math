from sympy import *
x, y, z = symbols('x y z')

expr = x**2 + x*y
srepr(expr)
expr = sin(x*y)/2 - x**2 + 1/y
srepr(expr)