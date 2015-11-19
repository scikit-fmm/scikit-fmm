from sympy import summation, symbols, diff, MatrixSymbol, Matrix
import sympy.printing.python as python
from sympy.abc import x,y,z,i,j,k,l
import re

a = symbols('a:64')

b0, b1, b2 = symbols('b:3')  # the known gridpoint
c0, c1, c2 = symbols('c:3')  # the unknown point on the zero level-set

expr = 0
l = 0
for i in range(4):
    for j in range(4):
        for k in range(4):
            expr += a[l] * x**i * y**j * z**k
            l +=1

print "interpolant expression"
print expr

dpdx_c = diff(expr,x).subs(x,c0).subs(y,c1).subs(z,c2)
dpdy_c = diff(expr,y).subs(x,c0).subs(y,c1).subs(z,c2)
dpdz_c = diff(expr,z).subs(x,c0).subs(y,c1).subs(z,c2)

eq0 = expr.subs(x,c0).subs(y,c1).subs(z,c2)

grad_p = Matrix((dpdx_c, dpdy_c, dpdz_c))
eq1 = grad_p.cross(Matrix((b0-c0,b1-c1,b2-c2)))

# eq1 = dpdx_c * (b1-c1) - dpdy_c * (b0-c0)

print
print
print re.sub(r"a(\d+)", r"a[\1-1]", str(eq1[0]))
print
print
print re.sub(r"a(\d+)", r"a[\1-1]", str(eq1[1]))
print
print
print re.sub(r"a(\d+)", r"a[\1-1]", str(eq1[2]))
print

# print "Jacobian elements"

# print diff(eq0, c0)
# print
# print diff(eq0, c1)

# print
# print diff(eq1, c0)
# print
# print diff(eq1, c1)
