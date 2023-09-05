import openturns as ot



def checkCovariance(sample, cov):
    cov_est = ot.NonStationaryCovarianceModelFactory().build(sample)
    MCov = cov.discretize(sample.getMesh())
    MCov_est = cov_est.discretize(sample.getMesh().getVertices())
    cov_L2 = MCov.computeGram().computeTrace()
    err_L2 = (MCov_est - MCov).computeGram().computeTrace()
    return (err_L2 / cov_L2)**0.5
    
    
