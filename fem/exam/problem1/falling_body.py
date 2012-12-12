from numpy import *
import matplotlib.pyplot as plt

g   = 9.81		# gravity
mu  = 1.8e-5	# viscosity in air
#mu  = 1.0
rho = 0.79		# density in air
B   = 1.0 		# global variable for stokes air resistance
A   = 1.0 		# gravity term 
volume = 1.0    # Volume
d   = 1.0       # global variable diameter

method = 0      # Forward Euler
#method = 1     # Backward Euler
#method = 2     # Crank-Nicolson

def read_command_line():
    parser = define_command_line_options()
    args = parser.parse_args()

    print 'T={}, dt={}, v0={}, method={}, makeplot={}'.format(
        args.T, args.dt, args.v0, args.method, args.makeplot)
    return args.T, args.dt, args.v0, args.method, args.makeplot

def define_command_line_options():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--T', '--stop_time', type=float,
                        default=30.0, help='end time of simulation [s]',
                        metavar='T')
    parser.add_argument('--dt', type=float,
                        default=1.0, help='Timestep [s]',
                        metavar='dt')
    parser.add_argument('--v0', type=float,
                        default=0.0, help='Initial velocity [m/s]',
                        metavar='v0')
    parser.add_argument('--method', type=int,
                        default=0, help='Solver method, 0: FE, 1: BE, 2: CN',
                        metavar='method')
    parser.add_argument('--makeplot', action='store_true',
                        help='Display plot [true/false]')
    return parser

def stepStokes(v,dt):
    if method == 0: return v*(1 + dt*B) + dt*A
    if method == 1: return (dt*A + v)/(1 - dt*B)
    if method == 2: return (A*dt + v*(1 + 0.5*B*dt))/(1 - 0.5*B*dt)

def step(v,dt):
	reynold = rho*d*abs(v)/mu

	return stepStokes(v,dt)

def solver(T, dt, v0, _method):
    '''
    Solves the vertical motion in a fluid from t=0 to t=T
    with timestep dt, diameter d, object density rho_body,
    fluid density rho, viscosity mu and initial velocity v0
    '''
    global A, B, d, rho, rho_body, volume, method
    
    method = _method
    d = 0.1
    rho_body = 1.79
    rho = 0.79
    volume = 1.0/24*pi*d**3

    A = g*(rho/rho_body - 1) # set global variable b based on object density
    B = -3*pi*d*mu/(rho_body*volume)
    #13.8117283951
    
    N = int(T/dt)+1
    v = zeros(N)
    v[0] = v0

    for n in range(1,N):
        v[n] = step(v[n-1],dt)
    return v

def main():
    T, dt, v0, method, makeplot = read_command_line()
    N = int(T/dt)+1
    t = linspace(0,T,N)
    colors = ['r','b','g']
    v_max = zeros(3)

    for method in [1,2]:
        v = solver(T, dt, v0, method)
        print "Simulated v_max, method %d=%g" % (method,v[-1])

        if makeplot:
            plt.plot(t,v,colors[method])
            plt.hold('on')

    v_e = volume*g*(rho - rho_body)/(3*pi*d*mu)
    print "Analytical v_max=",v_e
    if makeplot:
        v_e = A/B*(exp(B*t) - 1)
        plt.plot(t,v_e,'k')
        #plt.legend(["Forward Euler", "Backward Euler", "Crank-Nicolson", "Exact solution"])
        plt.legend(["Backward Euler", "Crank-Nicolson", "Exact solution"])

    plt.show()

if __name__ == '__main__':
    main()