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
    
    
def createMeshFromNumberOfIntervalWithBounds(N, dim, lowerbound, upperbound):
    temp = N**(1/dim)
    parameter = int(temp)
    
    
    discretization = [parameter] * dim
    mesher = ot.IntervalMesher(discretization)
    interval = ot.Interval(lowerbound, upperbound)
    
    return mesher.build(interval)
    

#création d'un maillage du pavé [0,1]^dim ayant à peu près (N+1) sommets
def createMeshFromNumberOfInterval(N, dim):
    lowerbound = [0]*dim
    upperbound = [1]*dim
    
    return createMeshFromNumberOfIntervalWithBounds(N, dim, lowerbound, upperbound)

def getTimeOfARealization(process):
    time1 = time.time()
    field = process.getRealization()
    time2 = time.time()
    
    return (time2 - time1), field 
    
    
def getMeanTimeAndErrorOfASample(myprocess, sampleSize):
    time1 = time.time()
    sample = myprocess.getSample(sampleSize)
    time2 = time.time()
    covModel = myprocess.getCovarianceModel()
    error = checkCovariance(sample, covModel)
    return ((time2 - time1)/sampleSize), error 
    
    
    
#création d'un graphique via matplotlib et le sauvegarde    
def create_figure_and_save(title,xlabel,ylabel,ydata,xdata,filename):
    fig, ax = plt.subplots(1,1)
    ax.bar(height=ydata, x=xdata) #ydata: coordonnées y | xdata :coordonnées x
    fig.suptitle(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    for i in range(len(ydata)):
        ax.text(i, ydata[i], round(ydata[i],5))
    
    plt.savefig(filename)
    