import openturns as ot
import openturns.viewer as otv
import matplotlib.pyplot as plt
import os
from spectralModels import *
from mySpectralGaussianProcess3D import *


ot.RandomGenerator.SetSeed(os.getpid())

fmax3 = [10.0]*3 
a = 1.0
smodel3 = SpectralModelofExponential3D1D(a) #modèle de covariance associé: t = [t1 t2 t3]  ->  exp[-a norm[t]]


dim = 3
size = 1000
N = int(size**(1/dim))
if N%2 == 1:
   N += 1
discretization = [N-1]*dim #discrétisation
lowerbound = [0.0]*dim
upperbound = [(N-1)/(2*fmax3[i]) for i in range(dim)] #longueur de chaque arête spatial dans chaque direction
mesher = ot.IntervalMesher(discretization)
interval = ot.Interval(lowerbound,upperbound)

mesh = mesher.build(interval)
process = MySpectralGaussianProcess3D(smodel3, mesh, discretization, upperbound)
field = process.getRealization()

field.exportToVTKFile("meshSpectral3D-1000.vtk")
