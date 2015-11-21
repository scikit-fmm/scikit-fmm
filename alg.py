from sympy import summation, symbols, diff, MatrixSymbol, Matrix
import sympy.printing.python as python
from sympy.abc import x,y,i,j

a = MatrixSymbol('a', 4, 4)

b0, b1 = symbols('b:2')  # the known gridpoint
c0, c1 = symbols('c:2')  # the unknown point on the zero level-set

expr = 0
for i in range(4):
    for j in range(4):
        expr += a[i,j] * x**i * y**j
print "interpolant expression"
print expr

dpdx_c = diff(expr,x).subs(x,c0).subs(y,c1)
dpdy_c = diff(expr,y).subs(x,c0).subs(y,c1)

eq0 = expr.subs(x,c0).subs(y,c1)

grad_p = Matrix((dpdx_c, dpdy_c, 0))
eq1 = grad_p.cross(Matrix((b0-c0,b1-c1,0)))[2]

# eq1 = dpdx_c * (b1-c1) - dpdy_c * (b0-c0)

print python(eq0)
print
print python(eq1)

# print "Jacobian elements"

# print diff(eq0, c0)
# print
# print diff(eq0, c1)

# print
# print diff(eq1, c0)
# print
# print diff(eq1, c1)