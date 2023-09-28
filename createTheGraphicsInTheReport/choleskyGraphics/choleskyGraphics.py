import openturns as ot
import openturns.viewer as otv
import matplotlib.pyplot as plt
import os

ot.RandomGenerator.SetSeed(os.getpid())


N = 100
dim = 1
mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))

g = mesh.draw()
v = otv.View(g)
plt.savefig("CholeskyDim1-100.jpg")

###########################################
dim=2
N= int(1000**(1/dim))
mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))

g = mesh.draw()
v = otv.View(g, (600,600), square_axes=True)
plt.savefig("CholeskyDim2-961.jpg")

###########################################
dim=3
N= int(1000**(1/dim))
mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))

mesh.exportToVTKFile("CholeskyMaillageDim3-729.vtk")


###########################################
N = 1001
dim = 1
mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))
myCovModel=ot.ExponentialModel([1.0], [1.0])
process = ot.GaussianProcess(myCovModel, mesh)
field = process.getRealization()

g = field.drawMarginal(0)
v=otv.View(g)
plt.savefig("ReaCholesky1D-1001nodes.jpg")




