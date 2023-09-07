import openturns as ot
import time
from testTool import checkCovariance

for dim in [1, 2, 3]:
    myCovModel = ot.ExponentialModel([1.0]*dim, [1.0]) #modèle de covariance
    for size in [10, 100, 1000, 10000]:
        N = int(size**(1/dim)) #nombre de noeuds
        mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim)) #maillage
        t0 = time.time()
        t1 = t0
        num = 0
        while t1 - t0 < 5:
            num += 1
            process = ot.GaussianProcess(myCovModel, mesh) #processus gaussien
            field = process.getRealization()#réalisation du processus en les noeuds du maillage
            t1 = time.time()
        print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "meanTime=", (t1 - t0) / num, "s")
 
print("\n")  
for num in [250, 500, 750, 1000]: #nb de réalisations
    for dim in [1, 2, 3]:
       myCovModel=ot.ExponentialModel([1.0]*dim, [1.0])
       for size in [10, 100, 1000]:
         N = int(size**(1/dim)) #nombre de noeuds
         mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim)) #maillage
         process = ot.GaussianProcess(myCovModel, mesh) #processus gaussien
         t0 = time.time()
         sample=process.getSample(num) #échantillon de num réalisations du processus en les noeuds du maillage
         t1 = time.time()
         err=checkCovariance(sample, myCovModel) #calcul erreur L2
         print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "nb realisations=", num, "meanTime=", (t1 - t0) / num, "s", "errL2=",err)
