from sympy import summation, symbols, diff, MatrixSymbol, Matrix
import sympy.printing.python as python
from sympy.abc import x,y,z,i,j,k,l
import re
import numpy as np
import pylab as plt
import itertools

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


offsets = np.array(((0,0,0), (1,0,0), (0,1,0),
                    (1,1,0), (0,0,1), (1,0,1),
                    (0,1,1), (1,1,1)))
#offsets = list(itertools.product((0,1), (0,1), (0,1)))

matrix = []

# ok these are the first 8 rows...
for xo,yo,zo in offsets:
    row = np.zeros(64)
    print (xo,yo,zo), expr.subs(x,xo).subs(y,yo).subs(z,zo)
    v = str(expr.subs(x,xo).subs(y,yo).subs(z,zo)).replace("a","")
    v = v.replace("+", ",")
    v = "("+v+",)"
    for index in eval(v):
        row[index] = 1
    print row
    matrix.append(tuple(row))


# now dfdx of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    aa = expr.diff(x).subs(x,xo).subs(y,yo).subs(z,zo)
    print (xo,yo,zo), aa
    v = str(aa)
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    #print v1
    #print v2
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))


# now dfdy of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    # working but how to transform into a matrix row?
    v = str(expr.diff(y).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    #print v1
    #print v2
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))

# now dfdz of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    # working but how to transform into a matrix row?
    v = str(expr.diff(z).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    #print v1
    #print v2
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))

# now d2fdxdy of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    v = str(expr.diff(x).diff(y).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))

# now d2fdxdz of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    v = str(expr.diff(x).diff(z).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))


# now d2fdydz of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    v = str(expr.diff(y).diff(z).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))

# now d3fdxdydz of this expression
for xo,yo,zo in offsets:
    row = np.zeros(64)
    v = str(expr.diff(x).diff(y).diff(z).subs(x,xo).subs(y,yo).subs(z,zo))
    v1 = re.sub(r"a(\d+)", r"1", v)
    v1 = v1.replace("+", ",")
    v1 = "("+v1+",)"
    v1 = eval(v1)
    v2 = re.sub(r"\S*a(\d+)", r"\1", v)
    v2 = v2.replace("+", ",")
    v2 = "("+v2+",)"
    v2 = eval(v2)
    for mult, index in zip(v1,v2):
        row[index] = mult
    print row
    matrix.append(tuple(row))

barr = np.array(matrix)

binv = np.linalg.inv(barr)
print binv
