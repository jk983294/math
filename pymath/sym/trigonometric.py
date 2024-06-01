import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

x, y, z = sympy.symbols('x y z')

cos(acos(x))
asin(1)

# simplify
trigsimp(sin(x)**2 + cos(x)**2)
trigsimp(sin(x)**4 - 2*cos(x)**2*sin(x)**2 + cos(x)**4)
trigsimp(sinh(x)/tanh(x))

# apply the sum or double angle identities
expand_trig(sin(x + y))
expand_trig(tan(2*x))