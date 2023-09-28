import openturns as ot
import openturns.viewer as otv
import matplotlib.pyplot as plt
import os

ot.RandomGenerator.SetSeed(os.getpid())

ot.ResourceMap.SetAsString("HMatrix-ClusteringAlgorithm", "median") #type de découpage hiérarchique de l'espace
ot.ResourceMap.SetAsString("HMatrix-CompressionMethod", "AcaPlus") #moyen de compression des blocs admissibles
ot.ResourceMap.SetAsScalar("HMatrix-AdmissibilityFactor", 100) #facteur d'admissibilité
ot.ResourceMap.SetAsUnsignedInteger("HMatrix-MaxLeafSize",250) #N_leaf, taille maximale des feuilles d'un cluster tree
ot.ResourceMap.SetAsScalar("HMatrix-AssemblyEpsilon", 1e-7) #epsilon d'assemblage 
ot.ResourceMap.SetAsScalar("HMatrix-RecompressionEpsilon", 1e-7) #epsilon de recompression


N = 1001 #nb de noeuds
mesh = ot.IntervalMesher([N-1]).build(ot.Interval([-1.0], [1.0]))
myCovModel=ot.ExponentialModel([1.0], [1.0])

process = ot.GaussianProcess(myCovModel, mesh)
process.setSamplingMethod(ot.GaussianProcess.HMAT)
field = process.getRealization()

g = field.drawMarginal(0)
v=otv.View(g)
plt.savefig("hmatrixRea.jpg")






