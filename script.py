#import sys
from os import system
import numpy as np

E  = 8
runs = 10
Z = 8
T = 1
Z_min = 2
Z_max = 8
ITER = 300
ITER_LOCAL = 20
slack = 1
c = 1

# energy statistics part
#for x in range (0, 6):
#    for y in range (0, 25):
#        system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E+x,runs,Z,T,Z_min,Z_max,ITER,ITER_LOCAL,slack,c))

# parameter dependance part
#for x in np.arange (1, 10, 1):
#    var = x
#    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(var,runs,Z,T,Z_min,Z_max,ITER,ITER_LOCAL,slack,c))

for x in range (5,100,5):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,var,Z,T,Z_min,Z_max,ITER,ITER_LOCAL,slack,c))

for x in np.arange (4,16,0.25):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,var,T,Z_min,Z_max,ITER,ITER_LOCAL,slack,c))

for x in np.arange (0.15,3,0.15):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,var,Z_min,Z_max,ITER,ITER_LOCAL,slack,c))

for x in range (1,7,1):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,var,Z_max,ITER,ITER_LOCAL,slack,c))

for x in range (3,15,1):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,Z_min,var,ITER,ITER_LOCAL,slack,c))

for x in range (50,1000,50):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,Z_min,Z_max,var,ITER_LOCAL,slack,c))

for x in range (5,100,5):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,Z_min,Z_max,ITER,var,slack,c))

for x in np.arange (0.5,1.5,0.1):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,Z_min,Z_max,ITER,ITER_LOCAL,var,c))

for x in np.arange (0.5,1.5,0.1):
    var = x
    system ("python lucy_comparisons.py {} {} {} {} {} {} {} {} {} {}".format(E,runs,Z,T,Z_min,Z_max,ITER,ITER_LOCAL,slack,var))

