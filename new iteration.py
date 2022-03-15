import assemblymethod as am
import layermethod as lm
import numpy as np
import matplotlib.pyplot as plt

numberoflayers = 2

layer0 = lm.layer(0,2,0.5,0.3)
layer1 = lm.layer(2,3,0.5,0.3)

layers = [(layer0), (layer1)]