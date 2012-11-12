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
	comparison_plot(f, u, om, filename='sin_powers_results/sin_powers_k='+str(i)+'.pdf')
	i+=1

taylor = kill_high_order_terms(f.series(n=10), 9, x) # n=10 means error should be of order x**10 or less

comparison_plot(f, taylor, Omega[len(Omega) - 1], filename='sin_powers_results/sin_powers_taylor.pdf')

# We can easily see that the taylor expansion is really bad. It uses the derivative at x=0 as weights
# so the approximation is probably better around x=0 than the least squares.