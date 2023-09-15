import openturns as ot
import time
import openturns.viewer as otv
import matplotlib.pyplot as plt
import math
from P1interpolationGaussianProcess import *
from testTool import *
import os

ot.RandomGenerator.SetSeed(os.getpid())



f1 = ot.SymbolicFunction(["x","y"], ["1-sqrt(x^2+y^2)"])
f2 = ot.SymbolicFunction(["x","y","z"], ["1-sqrt(x^2+y^2+z^2)"])
f = [f1,f2]
toEliminate = [ot.Greater(),ot.Greater()] #sens à donner: ot.Less->inférieur, ot.Greater-->supérieur

#modèle de covariance

##dimension 2
sigma = 1.0
sigma2 = sigma**2
b1 = 1.0
b2 = 4.0
covModel2 = ot.GeneralizedExponential([b1,b2],[sigma],2)

##dimension 3
theta=1.0
sigma=1.0
covModel3 = ot.GeneralizedExponential([1.0/(theta**(1/2))]*3, [sigma] ,2) #modèle de covariance: t1 t2 t3 -> (sigma**2) * exp(-(||t||**2)/theta)

#regroupement de données 
nbNodesBE= 750 #nombre de noeuds de la boite englobante permettant la création du maillage Ma
covModels = [covModel2,covModel3]



for dim in [2,3]:
   for size in [10, 100, 1000, 10000]:
      N = int(size**(1/dim))
      xMin = [-1.0]*dim
      xMax = [1.0]*dim       
      I = ot.Interval(xMin, xMax) #limites d'une boite englobante (BE) (pour créer le maillage Ma)
      N2 = [int(nbNodesBE**(1/dim))]*dim #discrétisation de la BE 
      levelSet = ot.LevelSet(f[dim-2], toEliminate[dim-2], 1-(0.5**(1/dim))) #levelSet indique que l'on retire les noeuds n de la BE
                                                                         #qui ont une valeur f[dim-2](n) qui est toEliminate[dim-2]
                                                                         #(i.e inférieure ou supérieure) à 1-(0.5**(1/dim))
      t0 = time.time()
      t1 = t0
      num = 0
      while t1 - t0 < 5:
        num += 1
        mesh = ot.LevelSetMesher(N2).build(levelSet, I) #maillage Ma qui est la BE discrétisée par N2 et qui a en moins les noeuds à retirer selon le critère levelSet
        boundingBox = ot.IntervalMesher([N-1]*dim).build(I) #maillage dont la boîte associée contient le maillage Ma
        bBprocess = ot.GaussianProcess(covModels[dim-2], boundingBox) #simulation par cholesky sur boundingBox
        process = P1interpolationGaussianProcess(mesh,bBprocess) #processus d'interpolation P1 pour simuler sur Ma
        field = process.getRealization() #réalisation sur Ma
        t1 = time.time()
      print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(),"method=", "P1Interpolation+Cholesky", "nb_nodes_boundingBox=", boundingBox.getVerticesNumber() ,"meanTime=", str((t1 - t0) / num)+"s")
     
print("\n")     
      
for dim in [2,3]:
   for size in [10, 100, 500, 1000]:
        for num in [250,1000,4000,10000]:
          xMin = [-1.0]*dim
          xMax = [1.0]*dim 
          boxLength = xMax[0] - xMin[0]
          I = ot.Interval(xMin, xMax) 
          N2 = [int(nbNodesBE**(1/dim))]*dim 
          levelSet = ot.LevelSet(f[dim-2], toEliminate[dim-2], 1-(0.5**(1/dim))) 
          mesh = ot.LevelSetMesher(N2).build(levelSet, I)
      
          N = int(size**(1/dim)) +1
          boundingBox = ot.IntervalMesher([N-1]*dim).build(I)
          bBprocess = ot.GaussianProcess(covModels[dim-2], boundingBox) #cholesky
          process = P1interpolationGaussianProcess(mesh,bBprocess)
          sample =  process.getSample(num)
          err = checkCovariance(sample, covModels[dim-2])
          h_bB = math.sqrt(dim) * (boxLength/(N-1)) #diamètre de tous les dim-simplexes composant le maillage boundingBox       
          print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(),"method=", "P1Interpolation+Cholesky","nb_nodes_boundingBox=", boundingBox.getVerticesNumber(),"h_bB=",  h_bB, "nb_realisations=", num, "errL2=",err)
      

#simulation par Cholesky directement sur le maillage de simulation Ma (mesh)
#afin de comparer l'erreur L2 de M réalisations sur Ma par différentes P1-interpolations 
#et l'erreur L2 de M réalisations sur Ma par Cholesky   
print("")
for dim in [2,3]:
   for num in [250,1000,4000,10000]:
      xMin = [-1.0]*dim
      xMax = [1.0]*dim       
      I = ot.Interval(xMin, xMax) 
      N2 = [int(nbNodesBE**(1/dim))]*dim 
      levelSet = ot.LevelSet(f[dim-2], toEliminate[dim-2], 1-(0.5**(1/dim))) 
      mesh = ot.LevelSetMesher(N2).build(levelSet, I)
      meshProcess = ot.GaussianProcess(covModels[dim-2], mesh)
      sample =  meshProcess.getSample(M)
      err = checkCovariance(sample, covModels[dim-2])
      print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(),"method=", "Cholesky","nb_realisations=",num, "errL2=",err)
