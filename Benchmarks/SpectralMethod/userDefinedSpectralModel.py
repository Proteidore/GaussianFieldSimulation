import openturns as ot
import math


def getIndex(k,discretization):
      lastPos = len(discretization)-1
      index = k[lastPos]
      for i in range(lastPos,0,-1):
        index = (index*discretization[i-1]) + k[i-1]
      return index 
      
#explication dans le cas de la dimension 2:
#cette classe construit un modèle spectral S
#sur un pavé [-fmax1,fmax1] x [-fmax2,fmax2]
#dont chaque arête [-fmax_i,fmax_i] (i allant de 1 à 2)
#a été subdivisée en N_i intervalles
#de longueur deltaF_i et dont on note les centres
#des intervalles f^i_0, ..., f^i_(N_i -1). Le
#modèle spectral est défini par morceaux et est constant
#sur les pavés de la forme 
#[f^1_j - (deltaF_1 /2),f^1_j + (deltaF_1 /2)[ x [f^2_k - (deltaF_1 /2),f^2_k + (deltaF_2 /2)[
#où j va de 0 à N_1 -1 et k va de 0 à N_2 -1 (donc les pavés sont identifiés par leur couple (j,k))
#Pour parcourir l'ensemble {(j,k), où j dans {0,...,N_1-1} et k dans {0,...,N_2-1}}
#on procède ainsi: on fixe k à 0, puis j parcourt
#par ordre croissant {0,...,N_1-1} et on obtient les (0,0), ..., (N_1-1,0)
#puis on fixe k à 1 puis j parcourt
#par ordre croissant {0,...,N_1-1} et et on obtient les (0,1), ..., (N_1-1,1)
#ainsi de suite tant que k est strictement inférieur à N_2
#Par ce parcours, on peut assigner aux (k,j) un indice index allant de 0 à (N_1 * N_2) -1 
#correspondant au moment où l'indice (k,j) est visité (temps 0, temps 1, ...) et à
#partir de (k,j), N1 et N2 on peut retrouver index en faisant getIndex([k,j],[N_1,N2]) 

class UserDefinedSpectralModel(ot.SpectralModel):

    def __init__(self, l, outputDimension, h : ot.HermitianMatrixCollection):
       self.inputDimension = len(l)
       self.outputDimension = outputDimension
       self.freqGrids = l #c'est une liste de ot.RegularGrid
       self.hermMatrix = h
       
       N = 1 
       for i in range(self.inputDimension):
          N = N * (self.freqGrids[i].getVertices()).getSize()
          
       assert N == self.hermMatrix.getSize()
       
       
       iDim = self.inputDimension
       self.starters = [ self.freqGrids[i].getStart() for i in range(iDim)]
       self.N = [ self.freqGrids[i].getN() for i in range(iDim)]
       self.enders = [ self.freqGrids[i].getValue(self.N[i]-1) for i in range(iDim)]
       self.steps = [ self.freqGrids[i].getStep() for i in range(iDim)]
       
       
       
       
    def __call__(self, f):
       
       iDim = self.inputDimension
       pos = [0] * iDim
       starters = self.starters
       N = self.N
       enders = self.enders
       steps = self.steps
       freqGrids = self.freqGrids
       
       
       
       for i in range(iDim):
          if f[i] < (starters[i] - (0.5*steps[i])) or f[i] > (enders[i] + (0.5*steps[i])):
             assert 1 == 0
          
          start, end = 0, N[i]-1 #dichotomie
          while(start != end):
             c = (start + end)//2
             if f[i] > (freqGrids[i].getValue(c) - (0.5*steps[i])) and f[i] <= (freqGrids[i].getValue(c) + (0.5*steps[i])):
                start, end = c, c
               
             else:
                if f[i] > (freqGrids[i].getValue(c) - (0.5*steps[i])):
                   start = c+1

                if f[i] <= (freqGrids[i].getValue(c) + (0.5*steps[i])):
                   end = c-1                   
                
          pos[i] = start  
          
          
       index = getIndex(pos,N)
       return self.hermMatrix.at(index)
         
    
    def getInputDimension(self):
          return self.inputDimension
          
    def getOutputDimension(self):
          return self.outputDimension
          
       
        
       
