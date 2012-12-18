import odespy
import numpy as np
import matplotlib.pyplot as plt
import sys

case = [[0.001,1.0], [0,1,3]]   # RK best, small time step
#case = [[0.4,15], [0,1,3]]   # BE best, large time step
#case = [[0.1,3], [0,1,2,3]]   # Runge kutta good, BE worst
#case = [[1.0,3], [0,1,3]]   # Runge kutta problems
#case = [[1.0,10], [0,1,3]]  # Runge kutta problems 2
#case = [[0.5,20], [1,2,3]]  # FE bye bye
#case = [[0.6,5], [1,3]]  # FE bye bye

solver_indices = case[1]
dt = case[0][0]
T = case[0][1]

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

dt = float(sys.argv[1]) if len(sys.argv) >= 2 else dt
T = float(sys.argv[2]) if len(sys.argv) >= 3 else T

N = int(round(T/dt))
t = np.linspace(0, N*dt, N+1)

# Compare with exact solution plotted on a very fine mesh
t_fine = np.linspace(0, T, 10001)
u_e = I*np.exp(B*t_fine) - A/B

solvers = [
    odespy.RK4(f),
    odespy.BackwardEuler(f, nonlinear_solver='Newton'),
    odespy.ForwardEuler(f, nonlinear_solver='Newton'),
    odespy.CrankNicolson(f, nonlinear_solver='Newton')]
legends = []
for solver_index in solver_indices:
    solver = solvers[solver_index]

    solver.set_initial_condition(0)
    u, t = solver.solve(t)
    plt.plot(t, u)
    plt.hold('on')
    legends.append(solver.__class__.__name__)
    if solver_index == 0: print "error_RK4=",abs(u[-1]-u_e[-1])
    if solver_index == 2: print "error_FE=",abs(u[-1]-u_e[-1])
    if solver_index == 1: print "error_BE=",abs(u[-1]-u_e[-1])
    if solver_index == 3: print "error_CN=",abs(u[-1]-u_e[-1])

plt.plot(t_fine, u_e, '-') # avoid markers by specifying line type

print "Stability criteria FE: ",(1 + B*dt)
print "Stability criteria CN: ",(1 + 0.5*B*dt)
print "dt-stability FE: ",-1./B
print "dt-stability FE (osc): ",-2./B
print "dt-stability CN: ",-2./B

legends.append('exact')
plt.legend(legends)
plt.title('Time step: %g' % dt)
plt.show()