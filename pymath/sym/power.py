import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

x, y = symbols('x y', positive=True)
a, b = symbols('a b', real=True)
z, t, c = symbols('z t c')

powsimp(x**a*x**b)
powsimp(x**a*y**a)
powsimp(t**c*z**c)
powsimp(t**c*z**c, force=True)
powsimp(z**2*t**2)
powsimp(sqrt(x)*sqrt(y))

expand_power_exp(x**(a + b))
expand_power_base((x*y)**a)
expand_power_base((z*t)**c)
powdenest((z**a)**b, force=True)