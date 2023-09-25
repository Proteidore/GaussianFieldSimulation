import openturns as ot
import time
from testTool import checkCovariance

ot.ResourceMap.SetAsString("HMatrix-ClusteringAlgorithm", "median") #type de découpage hiérarchique de l'espace
ot.ResourceMap.SetAsString("HMatrix-CompressionMethod", "AcaPlus") #moyen de compression des blocs admissibles
ot.ResourceMap.SetAsScalar("HMatrix-AdmissibilityFactor", 100) #facteur d'admissibilité
ot.ResourceMap.SetAsUnsignedInteger("HMatrix-MaxLeafSize",250) #N_leaf, taille maximale des feuilles d'un cluster tree
ot.ResourceMap.SetAsScalar("HMatrix-AssemblyEpsilon", 1e-7) #epsilon d'assemblage 
ot.ResourceMap.SetAsScalar("HMatrix-RecompressionEpsilon", 1e-7) #epsilon de recompression


for dim in [1, 2, 3]:
    for size in [10, 100, 1000, 10000]:
        N = int(size**(1/dim)) #nombre de noeuds
        mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim))#maillage
        t0 = time.time()
        t1 = t0
        num = 0
        while t1 - t0 < 5:
            num += 1
            process = ot.GaussianProcess(ot.ExponentialModel([1.0]*dim, [1.0]), mesh)#processus gaussien
            process.setSamplingMethod(ot.GaussianProcess.HMAT) #choix de simulation par H-matrices
            field = process.getRealization()  #réalisation du processus en les noeuds du maillage
            t1 = time.time()
        print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "meanTime=", (t1 - t0) / num, "s")


print("\n")  

for dim in [1, 2, 3]: #dimension de la maille
   for size in [10, 100, 1000]: 
       N = int(size**(1/dim)) #nombre de noeuds
       for num in [250, 500, 750, 1000]: #nombre de réalisations
         mesh = ot.IntervalMesher([N-1]*dim).build(ot.Interval([-10.0]*dim, [10.0]*dim)) #maillage
         myCovModel=ot.ExponentialModel([1.0]*dim, [1.0]) #modèle de covariance
         process = ot.GaussianProcess(myCovModel, mesh) #processus gaussien
         process.setSamplingMethod(ot.GaussianProcess.HMAT) #choix de simulation par H-matrices
         t0 = time.time()
         sample=process.getSample(num) #échantillon de num réalisations du processus en les noeuds du maillage
         t1 = time.time()
         err=checkCovariance(sample, myCovModel)#calcul erreur L2
         print("dim=", dim, "nb nodes=", mesh.getVerticesNumber(), "nb realisations=", num, "meanTime=", (t1 - t0) / num, "s", "errL2=",err)        
