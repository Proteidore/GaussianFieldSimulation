import math
import cmath
import openturns as ot
from userDefinedSpectralModel import *
from tool import *



class estimationSpectrale2D:
 
    def __init__(self,filteringWindows):
        self.filteringWindows = filteringWindows
        
        
    
    def __computeSpecialExpo(self,m,nbVertices): #m est un tuple ici de préférence (cas dim sup à 1)
        temp = 0
        for j in range(2):
            temp += ( m[j]*(1.0 - (1.0/nbVertices[j])) ) 
            
        return cmath.exp((1j)*cmath.pi*temp)
        
        
        
    def __evaluateFWindows(self, tab):
       value1 = self.filteringWindows(tab[0])
       value2 = self.filteringWindows(tab[1])
       return value1 * value2
    
    
    def buildFromSample(self,sample,nbVerticesEachDirection,upperbound):
       #Etape 1 : définir les grilles fréquentielles
       inputDim = 2

       discretization = []
       for i in range(inputDim):
          discretization.append(nbVerticesEachDirection[i]-1)
          
       spatialSteps = createSpatialSteps(discretization,upperbound)
       
       newUpperbound, fmax = [], []
       for i in range(inputDim):
          newUpperbound.append(upperbound[i]+spatialSteps[i])
          fmax.append(0.5/spatialSteps[i])
          
       spectralSteps= createSpectralSteps(newUpperbound)
       
       freqGrids=[]
       outputDim = sample.getDimension()
       
       for i in range(inputDim):
          freqGrids.append(ot.RegularGrid(-fmax[i]+(0.5*spectralSteps[i]),spectralSteps[i],nbVerticesEachDirection[i]))
          
       #Etape 2 : définir les matrices hermitiennes qui définiront par morceaux l'estimation de la densité spectrale
       
       N1=nbVerticesEachDirection[0]
       N2=nbVerticesEachDirection[1]
       N = N1 * N2
       deltaT = spatialSteps[0] * spatialSteps[1]
       cst = math.sqrt(deltaT/N)
       size = sample.getSize()
       fft_algorithm = ot.FFT()
     
       hermMatrixColl = []
       for i in range(N):
          hermMatrixColl.append(ot.HermitianMatrix(outputDim))
          
       for l in range(size):
           TFDMatrixColl = []
           
           for p in range(outputDim):
              matrix = ot.ComplexMatrix(N1,N2)
              pos = 0
              for j in range(N2):
                 for i in range(N1):
                    temp = cst * self.__evaluateFWindows((i/N1,j/N2)) *  self.__computeSpecialExpo((i,j),nbVerticesEachDirection)
                    matrix[i,j] = temp * sample[l][pos][p]
                    pos += 1
                    
              newMatrix = fft_algorithm.transform2D(matrix)
              TFDMatrixColl.append(newMatrix)
           
           pos = 0
           for j in range(N2):
               for i in range(N1):
                    vect = ot.ComplexMatrix(outputDim,1)
                    for p in range(outputDim):
                       vect[p,0] = TFDMatrixColl[p][i,j]
                    hermMatrixColl[pos] += ot.HermitianMatrix(vect*vect.conjugateTranspose())
                    pos += 1
         
                    
       for pos in range(N):
          hermMatrixColl[pos] = (1.0/size) * hermMatrixColl[pos]
       finalColl = ot.HermitianMatrixCollection(hermMatrixColl)
       
       return UserDefinedSpectralModel(freqGrids,outputDim, finalColl)
          
       
           
            
                   
       
