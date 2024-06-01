import sympy
from sympy import *

sympy.init_printing(use_unicode=True)

# compute derivatives, integrals, and limits, solve equations
x, y, t, z, nu = sympy.symbols('x y t z nu')

# derivative
sympy.diff(sin(x)*exp(x), x)
diff(exp(x**2), x)
diff(x**4, x, x, x)
diff(x**4, x, 3) # third derivative respect to x
diff(exp(x*y*z), x, 1, y, 2, z, 3)
exp(x*y*z).diff(x, 1, y, 2, z, 3)

# integrals
sympy.integrate(exp(x)*sin(x) + exp(x)*cos(x), x)
sympy.integrate(sin(x**2), (x, -oo, oo))
sympy.integrate(sin(x), (x, 0, pi / 2))
integrate(exp(-x), (x, 0, oo))
integrate(exp(-x**2 - y**2), (x, -oo, oo), (y, -oo, oo))
Integral(log(x)**2, x).doit()
Integral(sin(x**2), x).doit()

# limits
sympy.limit(sin(x)/x, x, 0)
Limit((cos(x) - 1)/x, x, 0).doit()
limit(1/x, x, 0, '+')
limit(1/x, x, 0, '-')


# Series Expansion
sympy.series(cos(x), x)
exp(sin(x)).series(x, 0, 4)  # expand at 0 with order 4 approximate
exp(x - 6).series(x, x0=6)

