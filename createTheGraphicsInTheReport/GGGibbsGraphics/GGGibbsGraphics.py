import openturns as ot
import openturns.viewer as otv
import matplotlib.pyplot as plt
import os

ot.RandomGenerator.SetSeed(os.getpid())

ot.ResourceMap.SetAsUnsignedInteger('GaussianProcess-GibbsMaximumIteration', 300)

N = 100 #nb de noeuds
mesh = ot.IntervalMesher([N-1]).build(ot.Interval([-1.0], [1.0]))
myCovModel=ot.ExponentialModel([1.0], [1.0])

process = ot.GaussianProcess(myCovModel, mesh)
process.setSamplingMethod(ot.GaussianProcess.GALLIGAOGIBBS)
field = process.getRealization()

g = field.drawMarginal(0)
v=otv.View(g)
plt.savefig("gibbsRea.jpg")
