from sympy import *
from approx1D import *
from fe_approx1D_v1 import *

x = symbols('x')
phi = basis(d=1)
nodes = [0,pi/2,pi]
elements = [[0,1],[1,2]]
f = sin(x)

A,b = assemble(nodes,elements,phi,f,symbolic=False)
print 'Final answer: '
#print A
#print b
c = A.LUsolve(b)

print c