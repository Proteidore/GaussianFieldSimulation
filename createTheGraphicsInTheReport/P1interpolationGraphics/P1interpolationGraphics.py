import openturns as ot
import time
import openturns.viewer as otv
import matplotlib.pyplot as plt
import math
from P1interpolationGaussianProcess import *
import os

ot.RandomGenerator.SetSeed(os.getpid())



dim=3
f = ot.SymbolicFunction(["x","y","z"], ["1-sqrt(x^2+y^2+z^2)"])
nbNodesBE= 750
xMin = [-1.0]*dim
xMax = [1.0]*dim       
I = ot.Interval(xMin, xMax) 
N2 = [int(nbNodesBE**(1/dim))]*dim 
levelSet = ot.LevelSet(f, ot.Greater(), 1-(0.5**(1/dim)))
mesh = ot.LevelSetMesher(N2).build(levelSet, I)

mesh.exportToVTKFile("meshP1Dim3-750.vtk")

####################################################
dim = 2
f = ot.SymbolicFunction(["x","y"], ["1-sqrt(x^2+y^2)"])


sigma = 1.0
sigma2 = sigma**2
b1 = 1.0
b2 = 4.0
covModel2 = ot.GeneralizedExponential([b1,b2],[sigma],2)


nbNodesBE= 750
size = 1000

N = int(size**(1/dim))
xMin = [-1.0]*dim
xMax = [1.0]*dim       
I = ot.Interval(xMin, xMax) 
N2 = [int(nbNodesBE**(1/dim))]*dim 
levelSet = ot.LevelSet(f, ot.Greater(), 1-(0.5**(1/dim)))
mesh = ot.LevelSetMesher(N2).build(levelSet, I)
boundingBox = ot.IntervalMesher([N-1]*dim).build(I) 
bBprocess = ot.GaussianProcess(covModel2, boundingBox)

process = P1interpolationGaussianProcess(mesh,bBprocess)
field = process.getRealization()

g = mesh.draw()
v = otv.View(g, (600,600), square_axes=True)


plt.savefig("meshP1Dim2-750.pdf")


field.exportToVTKFile("meshP1Dim2-750-961.vtk")

"""
g = field.draw()
g.setTitle(f"Réalisation par P1-interpolation")
view = otv.View(g, (600,600), square_axes=True)
plt.show()
"""
