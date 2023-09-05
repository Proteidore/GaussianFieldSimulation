import math
import cmath
import openturns as ot
from userDefinedSpectralModel import *



class estimationSpectrale1D:
 
    def __init__(self,filteringWindows):
        self.filteringWindows = filteringWindows 
    
    def __computeSpecialExpo(self,m,nbVertices): #m et nbVertices sont des entiers
        temp =  m*(1.0 - (1.0/nbVertices))  
        return cmath.exp((1j)*cmath.pi*temp)
        
        
        
    def __evaluateFWindows(self, value):
       return self.filteringWindows(value)
      
    
    
    def buildFromSample(self,sample,nbVertices,upperbound):
       #Etape 1 : définir les grilles fréquentielles
       inputDim = 1
       discretization = nbVertices - 1
       spatialStep = upperbound/discretization
       
       newUpperbound = upperbound + spatialStep 
       fmax = 0.5/spatialStep
       spectralStep = 1.0/newUpperbound
       outputDim = sample.getDimension()
       freqGrid = ot.RegularGrid(0.5*spectralStep,spectralStep,nbVertices//2)
          
       #Etape 2 : définir les matrices hermitiennes qui définiront par morceaux l'estimation de la densité spectrale
    
       N = nbVertices
       deltaT = spatialStep
       cst = math.sqrt(deltaT/N)
       size = sample.getSize()
       fft_algorithm = ot.FFT()
     
       hermMatrixColl = []
       for i in range(N//2):
          hermMatrixColl.append(ot.HermitianMatrix(outputDim))
       
       for l in range(size):
           TFDComplexColl = []
           
           for p in range(outputDim):
              complexColl = ot.ComplexCollection(N)
              for i in range(N):
                  temp = cst * self.__evaluateFWindows(i/N) *  self.__computeSpecialExpo(i,N)
                  complexColl[i] = temp * sample[l][i][p]
                    
              newComplexColl = fft_algorithm.transform(complexColl)
              TFDComplexColl.append(newComplexColl)
           
      
           for i in range(N//2):
              vect = ot.ComplexMatrix(outputDim,1)
              for p in range(outputDim):
                     vect[p,0] = TFDComplexColl[p][(N//2) + i]
                  
              hermMatrixColl[i] += ot.HermitianMatrix(vect*vect.conjugateTranspose())
         
                    
       for pos in range(N//2):
          hermMatrixColl[pos] = (1.0/size) * hermMatrixColl[pos]
       finalColl = ot.HermitianMatrixCollection(hermMatrixColl)
       
       return ot.UserDefinedSpectralModel(freqGrid, finalColl) 
  
