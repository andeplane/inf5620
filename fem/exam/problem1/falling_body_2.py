import odespy
import numpy as np
import matplotlib.pyplot as plt
import sys

g   = 9.81          # gravity
mu  = 1.0           # Viscosity
rho = 0.79          # Air density
rho_body = 1.79     # Object density
r   = 0.62035049089 # radius making volume = 1.0
volume = 4./3*np.pi*r**3

A = g*(rho/rho_body - 1) # set global variable b based on object density
B = -6*np.pi*r*mu/(rho_body*volume)

I = A/B

def f(u, t):
    return A + B*u

dt = float(sys.argv[1]) if len(sys.argv) >= 2 else 0.75
T = float(sys.argv[2]) if len(sys.argv) >= 3 else 10

N = int(round(T/dt))
t = np.linspace(0, N*dt, N+1)
solvers = [
    odespy.RK4(f),
    odespy.BackwardEuler(f, nonlinear_solver='Newton'),
    odespy.ForwardEuler(f, nonlinear_solver='Newton'),
    odespy.CrankNicolson(f, nonlinear_solver='Newton')]
legends = []
for solver in solvers:
    solver.set_initial_condition(0)
    u, t = solver.solve(t)
    plt.plot(t, u)
    plt.hold('on')
    legends.append(solver.__class__.__name__)
# Compare with exact solution plotted on a very fine mesh
t_fine = np.linspace(0, T, 10001)
u_e = I*np.exp(B*t_fine) - A/B
plt.plot(t_fine, u_e, '-') # avoid markers by specifying line type

print "Stability criteria FE: ",(1 + B*dt)
print "Stability criteria CN: ",(1 + 0.5*B*dt)

legends.append('exact')
plt.legend(legends)
plt.title('Time step: %g' % dt)
plt.show()