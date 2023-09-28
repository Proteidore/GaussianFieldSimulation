import cmath
import openturns as ot
import numpy as np
from tool import *

class MySpectralGaussianProcess3D:
    def __init__(self,spectralModel, spatialMesh, discretization, upperbound):
        dim = 3
        assert spatialMesh.getDimension() == dim
        assert spectralModel.getInputDimension() == dim 
        
        self.spectralModel = spectralModel
        self.spatialMesh = spatialMesh #c'est le maillage effectif (donc attention à discretization, upperbound et spatialSteps )
        self.inputDimension= spatialMesh.getDimension() 
        self.outputDimension = spectralModel.getOutputDimension()
        self.choleskyMatrices = []
        self.fft_algorithm = ot.FFT()
        
        self.upperbound=list(upperbound) 
        self.discretization=list(discretization) #là discretization désigne le nombre d'intervalles dans chaque direction de l'espace
        self.spatialSteps= createSpatialSteps(discretization,upperbound)
        self.fmax = [] #désigne la fréquence maximale dans chaque arête spectrale#désigne la fréquence maximale dans chaque arête spectrale
        
        for i in range(dim):
            assert self.discretization[i]%2 == 1
            self.discretization[i] += 1  #maintenant discretization désigne le nombre de points présents dans chaque direction de l'espace
            self.fmax.append(0.5/self.spatialSteps[i])
         
        
        for i in range(dim):
                self.upperbound[i] += self.spatialSteps[i] #il s'agit maintenant de l'upperbound du maillage non effectif associé à spatialMesh 
                
                
        self.nb_effective_vertices = self.discretization[0]*self.discretization[1]*self.discretization[2]
        self.spectralSteps= createSpectralSteps(self.upperbound)
        
        self.alpha = ot.ComplexCollection(self.nb_effective_vertices)
        self.alphaComputed = False




    def __computeSpecialExpo(self,m): #m est un tuple ici de préférence (cas dim sup à 1)
        temp = 0
        d = self.inputDimension
        for j in range(d):
            temp += ( m[j]*(1.0 - (1.0/self.discretization[j])) ) 
            
        return cmath.exp((-1j)*cmath.pi*temp)
        
        
    def __computeAlpha(self):
        if self.alphaComputed :
          return None
       
        self.alphaComputed = True
       
        deltaF = 1.0
        for i in range(len(self.spectralSteps)):
            deltaF *= self.spectralSteps[i]
             
        
        pos = 0
        factor = self.nb_effective_vertices * cmath.sqrt(deltaF)
        
        for k in range(self.discretization[2]):
           for j in range(self.discretization[1]):
              for i in range(self.discretization[0]):
                self.alpha[pos] = factor * self.__computeSpecialExpo((i,j,k))
                pos += 1
       
    def __getFk(self, k):
        l = []
        l.append((0.5+k[0])*self.spectralSteps[0])
        for j in range(1,len(k)):
           l.append( (0.5 + k[j])*self.spectralSteps[j] - self.fmax[j])
        
        return l

    def __createCholeskyMatrices(self):
        if len(self.choleskyMatrices) != 0 :
           return None
          
        for k in range(self.discretization[2]):
          for j in range(self.discretization[1]):
             for i in range(self.discretization[0]//2):
                 hermMatrix = self.spectralModel(self.__getFk((i,j,k)))
                 self.choleskyMatrices.append(hermMatrix.computeCholesky())
         
     
        
    def __computeZpk(self):
        n=self.outputDimension
        N= self.nb_effective_vertices
        result = []
        distribution = ot.Normal(4)
        
        for p in range(n):
           result.append(ot.ComplexTensor(self.discretization[0],self.discretization[1],self.discretization[2]))
           
        nHalves =  self.discretization[0]//2     
        pos = 0
        for k in range(self.discretization[2]):
          for j in range(self.discretization[1]):
            for i in range(nHalves):
            
               
               left = ot.ComplexCollection(n)
               right = ot.ComplexCollection(n)
               
               for p in range(n):
                  sample = distribution.getRealization()
                  left[p] = complex(sample[0],sample[1])
                  right[p] = complex(sample[2],sample[3])
               
               resultLeft = ot.ComplexCollection(self.choleskyMatrices[pos] * left)
               resultRight = ot.ComplexCollection(self.choleskyMatrices[pos] * right)
               
               for p in range(n):
                  (result[p])[nHalves + i,j,k] = resultRight[p] #ce choix de la conjugaison permet d'appliquer la FFT inverseTransform
                  (result[p])[nHalves-1-i,self.discretization[1]-1-j,self.discretization[2]-1-k] = resultLeft[p].conjugate()
               
               pos+=1
               
        return result
            
    def __computeXNp(self,zpk):
        n =  self.outputDimension
        N= self.nb_effective_vertices 
        sample = ot.Sample(N, n)

        for p in range(n):
           pseudoYNp = self.fft_algorithm.inverseTransform3D(zpk[p])
           pos = 0
           for k in range(self.discretization[2]):
             for j in range(self.discretization[1]):
               for i in range(self.discretization[0]):
                 temp = self.alpha[pos]*pseudoYNp[i,j,k]
                 sample[pos,p] = temp.real
                 pos += 1
        
        
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