import openturns as ot
import math


def getIdentity(dim):
   matrix = ot.HermitianMatrix(dim)
   for i in range(dim):
      matrix[i,i] = 1.0
   
   return matrix
    
def scalarProduct(v,w):
   n = 0
   for i in range(len(v)):
     n += (v[i]*w[i])
     
   return n  


#prend la valeur 1/2A sur [-A,A] et vaut 0 ailleurs
class UniformSpectralModel(ot.SpectralModel):
      
      def __init__(self, A : float):
         self.fmax = A
         self.constant = 0.5/self.fmax
           
      def __call__(self, f):
         matrix = ot.HermitianMatrix(1)
         if abs(f) <= self.fmax:
           matrix[0, 0] = self.constant
         else:
           matrix[0, 0] = 0
         
         return matrix
         
         
      def getInputDimension(self):
          return 1
          
      def getOutputDimension(self):
          return 1
         



# prend la valeur ((sigma**2) * (pi*theta)**(n/2)) * exp(-(pi**2) * theta * ||f||**2)*Identity(d)
# pour f dans R^inputDimension
# où n = inputDimension
# où d = outputDimension

#attention: la décroissance vers 0 à l'infini se veut très rapide
#par conséquent, plus le pas temporel se voudra petit, plus la plage de fréquence sera large, plus on risque de considérer l'image d'une grande fréquence
#comme la matrice nulle, ce qui peut être problématique
#il s'agit de la densité de la fonction d'autocovariance: t -> (sigma**2) * exp(-(||t||**2)/theta)
class GaussianSpectralModel(ot.SpectralModel):
      def __init__(self, inputDimension, outputDimension, theta, sigma):
         self.inputDimension = inputDimension 
         self.outputDimension = outputDimension
         self.theta = theta
         self.identity = getIdentity(self.outputDimension)
         temp = math.pi*theta
         self.value1 = (sigma**2) * (temp**(inputDimension/2.0))
         self.value2 = (- temp) * math.pi
         self.sigma = sigma


      def __call__(self, f): 
          norm2f2 = scalarProduct(f,f)
          
          temp = self.value2*norm2f2 #vérifier que les calculs sont bons en recalculant la TF
          temp = self.value1*math.exp(temp)
          
          return ((self.sigma**2) * temp) * self.identity
                
        
        
      def getInputDimension(self):
          return self.inputDimension
          
      def getOutputDimension(self):
          return self.outputDimension
 
class GaussianSpectralModelDim1(GaussianSpectralModel):
      def __init__(self, outputDimension, theta,sigma):
         super().__init__(1, outputDimension, theta,sigma)


      def __call__(self, f): 
          return super().__call__([f])          





#le modèle de covariance associé à cette densité spectrale qui va de R^3 dans R est
#t -> exp(-a*||t||) (où ||.|| est la norme euclidienne de R^3)       
class SpectralModelofExponential3D1D : 
       
       
       def __init__(self, a: float):
          self.inputDimension = 3 
          self.outputDimension = 1
          self.value = (2 * math.pi)**2
          self.identity = getIdentity(1)
       
          self.a2 = a**2
          self.cst = 8.0 * a * math.pi
       
       def __call__(self, f): 
          r2 =  scalarProduct(f,f)
          temp = self.a2 + (self.value * r2)
          temp = temp**2
          
          return (self.cst/temp)* self.identity
          
       def getInputDimension(self):
          return self.inputDimension
          
       def getOutputDimension(self):
          return self.outputDimension
          

#le modèle de covariance associé à cette densité spectrale qui va de R^2 dans R est
#(t1,t2) -> sigma2 * exp(-(t1/b1)^2 - (t2/b2)^2) 
class SpectralModelofAnExponential2D1D : 
       def __init__(self, sigma2, b1, b2):
          self.inputDimension = 2 
          self.outputDimension = 1
          self.identity = getIdentity(1)
          
          self.cons = sigma2 * b1 * b2 * math.pi 
          self.cons1 = (-1)* ((b1 * math.pi)**2)
          self.cons2 = (-1)* ((b2 * math.pi)**2)
       
         
       
       def __call__(self, f): 
          temp = math.exp((self.cons1 * (f[0]**2)) + (self.cons2 * (f[1]**2)))
          return (self.cons * temp) * self.identity
          
       def getInputDimension(self):
          return self.inputDimension
          
       def getOutputDimension(self):
          return self.outputDimension
 
