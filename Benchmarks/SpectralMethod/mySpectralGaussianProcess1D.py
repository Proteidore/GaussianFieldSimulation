import math
import openturns as ot
import numpy as np
from tool import *




#on va supposer que les spectralModel S en dimension 1 sont au format openturns normal. Donc si on fait 
#S(f), alors f est est un float
#et en dimension supérieure à 1, f est une liste de floats

class MySpectralGaussianProcess1D:
 
    def __init__(self,spectralModel, spatialMesh, discretization, upperbound):
        dim = 1
        assert spatialMesh.getDimension() == dim
        assert spectralModel.getInputDimension() == dim 
         
        
        self.spectralModel = spectralModel
        self.spatialMesh = spatialMesh #c'est le maillage effectif (donc attention à discretization, upperbound et spatialSteps )
        self.inputDimension= spatialMesh.getDimension() 
        self.outputDimension = spectralModel.getOutputDimension()
        self.choleskyMatrices = []
        self.fft_algorithm = ot.FFT()
        
        self.upperbound=list(upperbound) #il s'agit des temps maximaux sur chaque arête spatiale
        self.discretization=list(discretization) #il s'agit du nombre d'intervalles de même taille entre 0 et le temps maximal upperbound[0]
        self.spatialSteps= createSpatialSteps(discretization,upperbound)
        
        assert self.discretization[0]%2 == 1 #on impose que le nombre de points de simulation (donc qu'aussi le nombre d'intervalles de la partition fréquentielle)
                                             #soit paire
        self.discretization[0] += 1  #maintenant discretization désigne le nombre de points présents dans la direction x (de l'espace)
                                     #et désigne en même temps en combien d'intervalles l'arête spectrale a été subdivisée
                                     #Il y a un nombre paire d'intervalles ainsi                                    
        
         
        self.fmax = 0.5/self.spatialSteps[0] #désigne la fréquence maximale 
        self.upperbound[0] += self.spatialSteps[0] #il s'agit maintenant de l'upperbound du maillage non effectif associé à spatialMesh 
                
        self.nb_effective_vertices = self.discretization[0]
        self.spectralSteps= createSpectralSteps(self.upperbound) #pas de division du domaine fréquentiel
        
        self.alpha = ot.ComplexCollection(self.nb_effective_vertices)
        self.alphaComputed = False

    def __computeAlpha(self):
        if self.alphaComputed :
          return None
       
        self.alphaComputed = True
        deltaF = self.spectralSteps[0]
        N=self.nb_effective_vertices
        
        
        factor = N * math.sqrt(deltaF)
        for m in range(N):
          theta = math.pi*m*(1-(1/N))
          self.alpha[m] = factor * complex(math.cos(theta),-math.sin(theta))
       
    def __getFk(self, k):
        return (0.5+k)*self.spectralSteps[0]

    def __createCholeskyMatrices(self):
        if len(self.choleskyMatrices) != 0 :
          return None
          
        for k in range(self.nb_effective_vertices//2):
           hermMatrix = self.spectralModel(self.__getFk(k))
           self.choleskyMatrices.append(hermMatrix.computeCholesky()) 
      
    def __computeZpk(self):
        n=self.outputDimension
        N= self.nb_effective_vertices
        result = []
        distribution = ot.Normal(4)
        
        for p in range(n):
           result.append(ot.ComplexCollection(N))
         
        nHalves =  self.discretization[0]//2 
        for i in range(nHalves):
            
          left = ot.ComplexCollection(n)
          right = ot.ComplexCollection(n)
          for p in range(n):
            sample = distribution.getRealization()
            left[p] = complex(sample[0],sample[1])
            right[p] = complex(sample[2],sample[3])
               
          resultLeft = ot.ComplexCollection(self.choleskyMatrices[i] * left)
          resultRight = ot.ComplexCollection(self.choleskyMatrices[i] * right)
               
          for p in range(n):
            (result[p])[nHalves + i] = resultRight[p] #ce choix de la conjugaison permet d'appliquer la FFT inverseTransform
            (result[p])[nHalves-1-i] = resultLeft[p].conjugate()

               
        return result
    


     
    def __computeXNp(self,zpk):
        n =  self.outputDimension
        N= self.nb_effective_vertices 
        sample = ot.Sample(N, n)

        for p in range(n):
           pseudoYNp = self.fft_algorithm.inverseTransform(zpk[p])
           for i in range(N):
              temp = self.alpha[i]*pseudoYNp[i] #pseudo car il faut multiplier par N pour avoir YNp et self.alpha[i] contient cette multiplication
              sample[i,p] = temp.real
        
        
        return sample

    def getRealization(self):
        
        self.__computeAlpha()
        self.__createCholeskyMatrices()
        zpk = self.__computeZpk()
        xNp = self.__computeXNp(zpk)
        
        
        return ot.Field(self.spatialMesh,xNp) 

        
       
 
    def getSample(self, size):
        sample = ot.ProcessSample(self.spatialMesh,0, self.outputDimension)
        for i in range(size):
            sample.add(self.getRealization())
        return sample
        
    def getOutputDimension(self):
       return self.outputDimension
       
    def getMesh(self):
       return self.spatialMesh