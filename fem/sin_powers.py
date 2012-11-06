from approx1D import *
from fe_approx1D_v1 import *
from sympy import *
N = 4

x = Symbol('x')
k = Symbol('k')
f = sin(x)
phi = [x**(2*i+1) for i in range(N+1)]
Omega = [[0,k*pi/2] for k in range(2,13)]
#Omega = [[0,pi/2]]
i=2
for om in Omega:
	u = least_squares(f, phi, om)
	comparison_plot(f, u, om, filename='sin_powers_k='+str(i)+'.pdf')
	i+=1

# The taylor series of degree 9 is really bad. 