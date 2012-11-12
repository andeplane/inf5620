from sympy import *
from approx1D import *
from fe_approx1D_v1 import *

Omega = [0,2]
nodes = [0,1,1.2,2]
elements = [[0,1], [1,2], [2,3]]

nodes.insert(3,1.6)
elements.append([3,4])
print nodes
print elements
