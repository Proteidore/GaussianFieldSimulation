import math
import openturns as ot
import openturns.viewer as otv
import matplotlib.pyplot as plt
import time
from spectralModels import *
from mySpectralGaussianProcess1D import *
from mySpectralGaussianProcess2D import *
from mySpectralGaussianProcess3D import *
from testTool import *
from estimationSpectraleDim1 import *
from estimationSpectraleDim2 import *
from estimationSpectraleDim3 import *
import os

ot.RandomGenerator.SetSeed(os.getpid())


#création des modèles spectraux
##dimension 1
 
fmax1 = [5.0]
outputDim = 1 #dimension de sortie
smodel1= UniformSpectralModel(fmax1[0])


##dimension 2 
weight = 2.0
k1=40.00/weight #pulsation maximale de l'arete spectrale sur l'axe x
k2=10.0/weight #pulsation maximale de l'arete spectrale sur l'axe y
fmax2 = [k1/(2*math.pi),k2/(2*math.pi)]
sigma = 1.0
sigma2 = sigma**2
b1 = 1.0
b2 = 4.0
smodel2 = SpectralModelofAnExponential2D1D(sigma2, b1, b2) #modèle de covariance associé: t1 t2 -> sigma2 exp [ - [t1/b1]**2 - [t2/b2]**2  ]

##dimension 3
fmax3 = [10.0]*3 
a = 1.0
smodel3 = SpectralModelofExponential3D1D(a) #modèle de covariance associé: t = [t1 t2 t3]  ->  exp[-a norm[t]] 


#regroupement de données
fmax = [fmax1,fmax2,fmax3]
smodels = [smodel1,smodel2,smodel3]
spectralClasses= [MySpectralGaussianProcess1D,MySpectralGaussianProcess2D,MySpectralGaussianProcess3D]
estimateurs = [estimationSpectrale1D(ot.Hamming()),estimationSpectrale2D(ot.Hamming()),estimationSpectrale3D(ot.Hamming())]
compare = [compareSpectralModelsNormInfDim1,compareSpectralModelsNormInfDim2,compareSpectralModelsNormInfDim3] #liste de fonctions



for dim in [1,2,3]:
  for size in [10, 100, 1000, 10000]:
      N = int(size**(1/dim))
      if N%2 == 1:
        N += 1
      discretization = [N-1]*dim #discrétisation
      lowerbound = [0.0]*dim
      upperbound = [(N-1)/(2*fmax[dim-1][i]) for i in range(dim)] #longueur de chaque arête spatial dans chaque direction
      mesher = ot.IntervalMesher(discretization)
      interval = ot.Interval(lowerbound,upperbound)
      mesh = mesher.build(interval)
      t0 = time.time()
      t1 = t0
      num = 0
      while t1 - t0 < 5:
        num += 1
        process = spectralClasses[dim-1](smodels[dim-1], mesh, discretization, upperbound) #processus gaussien
        field = process.getRealization()#réalisation du processus en les noeuds du maillage
        t1 = time.time()
      print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "meanTime=", (t1 - t0) / num, "s")

print("\n")  

for dim in [1,2,3]:
   for num in [100,500,1000,2000]: #nb de réalisations
      for size in [100,500,1000]:
         N = int(size**(1/dim))
         if N%2 == 1:
           N += 1
         discretization = [N-1]*dim #discrétisation
         lowerbound = [0.0]*dim
         upperbound = [(N-1)/(2*fmax[dim-1][i]) for i in range(dim)]
         mesher = ot.IntervalMesher(discretization)
         interval = ot.Interval(lowerbound,upperbound)
         mesh = mesher.build(interval)
         spectralSteps = [(2*fmax[dim-1][i])/N for i in range(dim)]
         process = spectralClasses[dim-1](smodels[dim-1], mesh, discretization, upperbound) #processus gaussien
         
         sample = process.getSample(num)
         err = None
         if dim == 1:
            estimated_SM = estimateurs[dim-1].buildFromSample(sample,N,upperbound[0])
            err = compare[0](estimated_SM,smodels[0], -fmax[0][0] + 0.5*spectralSteps[0], spectralSteps[0], N)
         else:
            estimated_SM = estimateurs[dim-1].buildFromSample(sample,[N]*dim,upperbound)
            err = compare[dim-1](estimated_SM,smodels[dim-1], [-fmax[dim-1][i] + 0.5*spectralSteps[i] for i in range(dim)], spectralSteps, [N]*dim)
         print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "nb realisations=", num, "errSpecNormInf=",err)
         
         
         