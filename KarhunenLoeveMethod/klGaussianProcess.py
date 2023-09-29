import openturns as ot
import math
import time

class KLGaussianProcess:

     def __init__(self,mesh,algo,nbModes=0):  #attention!!!: covModel doit vérifier des propriétés

         self.mesh = mesh
         self.algo = algo
         self.outputDim = self.algo.getCovarianceModel().getOutputDimension()
         
         if nbModes > 0:
             self.algo.setNbModes(nbModes)   #remarque: le nb de eigenvalues calculé sera le min entre nbmodes et l'entier calculé à partir de threshold
             #sauf si threshold==0. Dans ce cas, nbmodes sera d'office sélectionné
             
         self.actualProcess = None
         self.nbEigenvalues = 0
         self.result = None


     def getRealization(self):
         if self.actualProcess == None:
             self.algo.run() #ce qui dure le plus longtemps quand le maillage et la base fonctionnelle deviennent élevés (dépend aussi du seuil s)
             self.result = self.algo.getResult()
             
             eigenvalues = self.result.getEigenvalues()
             for i in range(len(eigenvalues)):
                 if eigenvalues[i] == 0:
                     break
                 self.nbEigenvalues += 1
             self.actualProcess = ot.FunctionalBasisProcess(ot.Normal(self.nbEigenvalues), [self.result.getScaledModes()[i] for i in range(self.nbEigenvalues)], self.mesh)

         return self.actualProcess.getRealization()


     def getSample(self, size):
         if self.actualProcess == None:
             self.algo.run() #ce qui dure le plus longtemps quand le maillage et la base fonctionnelle deviennent élevés (dépend aussi du seuil s)
             self.result = self.algo.getResult()
             
             eigenvalues = self.result.getEigenvalues()
             for i in range(len(eigenvalues)):
                 if eigenvalues[i] == 0:
                     break
                 self.nbEigenvalues += 1
             self.actualProcess = ot.FunctionalBasisProcess(ot.Normal(self.nbEigenvalues), [self.result.getScaledModes()[i] for i in range(self.nbEigenvalues)], self.mesh)

         return self.actualProcess.getSample(size)

     def getCovarianceModel(self):
         return self.actualProcess.getCovarianceModel()

     def getNbOfModes(self):
         return self.nbEigenvalues

     def getMesh(self):
         return self.mesh

class KLP1AGaussianProcess(KLGaussianProcess):

   def  __init__(self,mesh, covModel, threshold, nbModes=0):
       algo = ot.KarhunenLoeveP1Algorithm(mesh, covModel, threshold)
       super().__init__(mesh,algo,nbModes)


class KLQAGaussianProcess(KLGaussianProcess):

   def  __init__(self, meshDomain, bounds, covariance, experiment, funcBasis, mustScale, s, nbModes=0):
       mesh = meshDomain.getMesh()
       xMin = mesh.getVertices().getMin()
       xMax = mesh.getVertices().getMax()
       algo = ot.KarhunenLoeveQuadratureAlgorithm(ot.Interval(xMin, xMax), bounds, covariance, experiment, funcBasis, mustScale, s)
       super().__init__(mesh,algo,nbModes)
