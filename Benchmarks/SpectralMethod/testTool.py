import openturns as ot
import time
import matplotlib.pyplot as plt


def checkCovariance(sample, cov):
    cov_est = ot.NonStationaryCovarianceModelFactory().build(sample)
    MCov = cov.discretize(sample.getMesh())
    MCov_est = cov_est.discretize(sample.getMesh().getVertices())
    cov_L2 = MCov.computeGram().computeTrace()
    err_L2 = (MCov_est - MCov).computeGram().computeTrace()
    return (err_L2 / cov_L2)**0.5
    
def createStationaryCovarianceModel(mesh, Cstat):
    covariance = ot.CovarianceMatrix(mesh.getVerticesNumber())
    for k in range(mesh.getVerticesNumber()):
        t = mesh.getVertices()[k]
        for ll in range(k + 1):
           s = mesh.getVertices()[ll]
           covariance[k, ll] = Cstat(t-s)
            
    return ot.UserDefinedCovarianceModel(mesh, covariance)
    
  
    
def getNormInfDim1(smodel1, freqGridStarter, spectralStep, N):
    assert smodel1.getInputDimension() == 1 
    
    maxValue = 0
    
    for i in range(N):
       freq = freqGridStarter + i*spectralStep
       temp = smodel1(freq) 
       value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
       if maxValue < value : 
          maxValue = value
          
    return maxValue**0.5

def compareSpectralModelsNormInfDim1(smodel1,smodel2, freqGridStarter, spectralStep, N):
    assert smodel1.getInputDimension() ==  smodel2.getInputDimension() == 1 
    assert smodel1.getOutputDimension() ==  smodel2.getOutputDimension()
    
    maxValue = 0
    
    for i in range(N):
       freq = freqGridStarter + i*spectralStep
       temp = smodel2(freq) - smodel1(freq) #calcul de la norme de Frobenius de la différence
       value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
       if maxValue < value : 
          maxValue = value
          
    return maxValue**0.5
    
    
def getNormInfDim2(smodel1, freqGridStarters, spectralSteps, N):
    assert smodel1.getInputDimension()  == 2 
    maxValue = 0
    for j in range(N[1]):
       freq2 = freqGridStarters[1] + j*spectralSteps[1]
       for i in range(N[0]):
          freq1 = freqGridStarters[0] + i*spectralSteps[0]
          temp = smodel1((freq1,freq2)) 
          value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
          if maxValue < value : 
             maxValue = value
          
    return maxValue**0.5     
    
    
def compareSpectralModelsNormInfDim2(smodel1,smodel2, freqGridStarters, spectralSteps, N):
    assert smodel1.getInputDimension() ==  smodel2.getInputDimension() == 2 
    assert smodel1.getOutputDimension() ==  smodel2.getOutputDimension()
    
    maxValue = 0
    
    for j in range(N[1]):
       freq2 = freqGridStarters[1] + j*spectralSteps[1]
       for i in range(N[0]):
          freq1 = freqGridStarters[0] + i*spectralSteps[0]
          temp = smodel2((freq1,freq2)) - smodel1((freq1,freq2)) #calcul de la norme de Frobenius de la différence
          value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
          if maxValue < value : 
             maxValue = value
          
    return maxValue**0.5    
    
    
def getNormInfDim3(smodel1, freqGridStarters, spectralSteps, N):
    assert smodel1.getInputDimension()  == 3 
    maxValue = 0
    for k in range(N[2]):
       freq3 = freqGridStarters[2] + k*spectralSteps[2]
       for j in range(N[1]):
          freq2 = freqGridStarters[1] + j*spectralSteps[1]
          for i in range(N[0]):
            freq1 = freqGridStarters[0] + i*spectralSteps[0]
            temp = smodel1((freq1,freq2,freq3)) 
            value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
            if maxValue < value : 
               maxValue = value
          
    return maxValue**0.5      
    
def compareSpectralModelsNormInfDim3(smodel1,smodel2, freqGridStarters, spectralSteps, N):
    assert smodel1.getInputDimension() ==  smodel2.getInputDimension() == 3 
    assert smodel1.getOutputDimension() ==  smodel2.getOutputDimension()
    
    maxValue = 0
    
    for k in range(N[2]):
       freq3 = freqGridStarters[2] + k*spectralSteps[2]
       for j in range(N[1]):
          freq2 = freqGridStarters[1] + j*spectralSteps[1]
          for i in range(N[0]):
            freq1 = freqGridStarters[0] + i*spectralSteps[0]
            temp = smodel2((freq1,freq2,freq3)) - smodel1((freq1,freq2,freq3)) #calcul de la norme de Frobenius de la différence
            value = ot.SquareMatrix((temp * temp).real()).computeTrace() #car hermitien
            if maxValue < value : 
               maxValue = value
          
    return maxValue**0.5  
 
#tiré d'openturns 
def graphOfTwoModelsInputDim1OutputDim1(smodel1,smodel2, freqGridStarter, spectralStep, N):
    assert smodel1.getInputDimension() ==  smodel2.getInputDimension() == 1 
    assert smodel1.getOutputDimension() ==  smodel2.getOutputDimension() == 1
    plotSample = ot.Sample(N, 3)

    
    for k in range(N):
       freq = freqGridStarter + k * spectralStep
       plotSample[k, 0] = freq
       plotSample[k, 1] = abs(smodel1(freq)[0, 0]) #modèle estimé
       plotSample[k, 2] = abs(smodel2(freq)[0, 0]) #modèle théorique

    
    graph = ot.Graph(
        "Estimated spectral function - Validation",
       "Frequency",
       "Spectral density function",
       True,
       "topright",
       1.0,
       ot.GraphImplementation.LOGY,
    )

    # The first curve is the estimate density as function of frequency
    curve1 = ot.Curve(plotSample.getMarginal([0, 1]))
    curve1.setColor("blue")
    curve1.setLegend("Estimate model")

    # The second curve is the theoritical density as function of frequency
    curve2 = ot.Curve(plotSample.getMarginal([0, 2]))
    curve2.setColor("red")
    curve2.setLegend("Theoretical model")

    graph.add(curve1)
    graph.add(curve2)
    return graph
    
