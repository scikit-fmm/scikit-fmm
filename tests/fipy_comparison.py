import fipy as fp
import numpy as np
from fipy.tools import serial
import skfmm


mesh = fp.Grid1D(dx = .5, nx = 8, communicator=serial)
value = (-1, -1, -1, -1, 1, 1, 1, 1)
var = fp.DistanceVariable(mesh = mesh, value=value)
var.calcDistanceFunction()
answer = (-1.75, -1.25, -.75, -0.25, 0.25, 0.75, 1.25, 1.75)
assert var.allclose(answer)
var2 = skfmm.distance(value, 0.5)
np.testing.assert_allclose(var2, answer)

dx = 1e-10
mesh = fp.Grid1D(dx = dx, nx = 8, communicator=serial)
var = fp.DistanceVariable(mesh = mesh, value=value)
var.calcDistanceFunction()
answer = np.arange(8) * dx - 3.5 * dx
assert var.allclose(answer)
var2 = skfmm.distance(value, dx)
np.testing.assert_allclose(var2, answer)

dx = 1.
dy = 2.
mesh = fp.Grid2D(dx = dx, dy = dy, nx = 2, ny = 3)
value = (-1, 1, 1, 1, -1, 1)
var = fp.DistanceVariable(mesh = mesh, value=value)

var.calcDistanceFunction()
vbl = -dx * dy / np.sqrt(dx**2 + dy**2) / 2.
vbr = dx / 2
vml = dy / 2.
crossProd = dx * dy
dsq = dx**2 + dy**2
top = vbr * dx**2 + vml * dy**2
sqrt = crossProd**2 *(dsq - (vbr - vml)**2)
sqrt = np.sqrt(max(sqrt, 0))
vmr = (top + sqrt) / dsq
answer = (vbl, vbr, vml, vmr, vbl, vbr)
assert var.allclose(answer)

# are the skfmm x and y backwards?
phi=np.array(value).reshape(3,2)
var2 = skfmm.distance(phi, [2.0, 1.0])
np.testing.assert_allclose(var2.flatten(), var.value)


#########
def fipy_test_extension_var():
    mesh = fp.Grid2D(dx = 1., dy = 1., nx = 2, ny = 2)
    var = fp.DistanceVariable(mesh = mesh, value = (-1, 1, 1, 1))
    var.calcDistanceFunction()
    extensionVar = fp.CellVariable(mesh = mesh, value = (-1, .5, 2, -1))
    tmp = 1 / np.sqrt(2)
    print var.value
    assert var.allclose((-tmp / 2, 0.5, 0.5, 0.5 + tmp))

var.extendVariable(extensionVar)
assert extensionVar.allclose((1.25, .5, 2, 1.25))
print extensionVar.value


mesh = fp.Grid2D(dx = 1., dy = 1., nx = 3, ny = 3)
var = fp.DistanceVariable(mesh = mesh, value = (-1, 1, 1,
                                               1, 1, 1,
                                               1, 1, 1))
var.calcDistanceFunction()
extensionVar = fp.CellVariable(mesh = mesh, value = (-1, .5, -1,
                                                   2, -1, -1,
                                                  -1, -1, -1))


v1 = 0.5 + tmp
v2 = 1.5
tmp1 = (v1 + v2) / 2 + np.sqrt(2. - (v1 - v2)**2) / 2
tmp2 = tmp1 + 1 / np.sqrt(2)
assert var.allclose((-tmp / 2, 0.5, 1.5, 0.5, 0.5 + tmp,
                      tmp1, 1.5, tmp1, tmp2))

answer = (1.25, .5, .5, 2, 1.25, 0.9544, 2, 1.5456, 1.25)
var.extendVariable(extensionVar)
assert extensionVar.allclose(answer, rtol = 1e-4)





#########

mesh = fp.Grid1D(dx = 1., nx = 3)
var = fp.DistanceVariable(mesh = mesh, value = (-1, 1, -1))
var.calcDistanceFunction()
assert var.allclose((-0.5, 0.5, -0.5))
np.testing.assert_allclose(var.value, skfmm.distance([-1, 1, -1]))

