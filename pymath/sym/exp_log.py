import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

x, y = symbols('x y', positive=True)
n = symbols('n', real=True)
z, t, c = symbols('z t c')

expand_log(log(x*y))
expand_log(log(x/y))
expand_log(log(x**2))
expand_log(log(x**n))
expand_log(log(z*t))
expand_log(log(z**2))
expand_log(log(z**2), force=True) # force can be used to ignore assumptions


logcombine(log(x) + log(y))
logcombine(n*log(x))
logcombine(n*log(z))
logcombine(n*log(z), force=True)