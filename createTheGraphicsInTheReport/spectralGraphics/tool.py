def createSpatialSteps(discretization,upperbound):
    l = []
    for i in range(len(upperbound)):
        v = upperbound[i]/discretization[i]
        l.append(v)
    return l

def createSpectralSteps(upperbound):
    l = []
    for i in range(len(upperbound)):
        v = 1.0/upperbound[i]
        l.append(v)
    return l

