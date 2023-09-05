import openturns as ot
import math


class P1interpolationGaussianProcess : 
      
      def __init__(self, mesh, envelopingProcess):
       self.mesh = mesh
       self.envelopingProcess = envelopingProcess
       self.outputDim = envelopingProcess.getOutputDimension()
       self.p1Interpolation = ot.P1LagrangeInterpolation(envelopingProcess.getMesh(),mesh,self.outputDim)
       
      def getRealization(self):
        preField = self.envelopingProcess.getRealization()
        return ot.Field(self.mesh, self.p1Interpolation.__call__(preField))
        
      def getSample(self, size):
        sample = ot.ProcessSample(self.mesh, 0, self.outputDim)
        for i in range(size):
           sample.add(self.getRealization())
           
        return sample
        
      def changeEnvelopingProcess(self, newEnvelopingProcess):
         self.envelopingProcess = newEnvelopingProcess
         self.outputDim = self.envelopingProcess.getOutputDimension()
         self.p1Interpolation = ot.P1LagrangeInterpolation(self.envelopingProcess.getMesh(),self.mesh,self.outputDim)
         

      def getMesh(self):
         return self.mesh
         
      def getOutputDimension(self):
         return self.outputDim
       
      def getEnvelopingMesh(self):
         return self.envelopingProcess.getMesh()
         
 
