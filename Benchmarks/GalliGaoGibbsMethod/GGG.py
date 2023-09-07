import openturns as ot
import time
from testTool import checkCovariance

ot.ResourceMap.SetAsUnsignedInteger('GaussianProcess-GibbsMaximumIteration', 300)

dim=1
myCovModel = ot.ExponentialModel([1.0]*dim, [1.0]) #modèle de covariance

for nbIt in [100,200,300]: #nombre d'itérations pour la méthode de Galli-Gao-Gibbs 
  ot.ResourceMap.SetAsUnsignedInteger('GaussianProcess-GibbsMaximumIteration', nbIt)
  for size in [10, 100]: #nombre de noeuds
     N = size
    
     t0 = time.time()
     t1 = t0
     num = 0
     while t1 - t0 < 5:
       num += 1
       mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim)) #maillage
       process = ot.GaussianProcess(myCovModel, mesh) #processus gaussien
       process.setSamplingMethod(ot.GaussianProcess.GALLIGAOGIBBS) #choix de simulation par Galli-Gao-Gibbs
       field = process.getRealization()
       t1 = time.time()
     print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(),"nb Iterations=", nbIt, "meanTime=", (t1 - t0) / num, "s")
 
print("\n")  
for nbIt in [10,50,100]: #nombre d'itérations pour la méthode de Galli-Gao-Gibbs 
  ot.ResourceMap.SetAsUnsignedInteger('GaussianProcess-GibbsMaximumIteration', nbIt)
  for num in [250, 500, 750, 1000]: #nombre de réalisations
     for size in [10, 100]: #nombre de noeuds
          N = size
          mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))
          process = ot.GaussianProcess(myCovModel, mesh) #processus gaussien
          process.setSamplingMethod(ot.GaussianProcess.GALLIGAOGIBBS) #choix de simulation par Galli-Gao-Gibbs
          t0 = time.time()
          sample=process.getSample(num) #échantillon de num réalisations du processus en les noeuds du maillage
          t1 = time.time()
          err=checkCovariance(sample, myCovModel) #calcul erreur L2
          print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "nb realisations=", num, "nb Iterations=", nbIt,"meanTime=", (t1 - t0) / num, "s", "errL2=",err)
