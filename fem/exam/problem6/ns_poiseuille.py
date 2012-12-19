from dolfin import *
import numpy as np

# Set parameter values
rho = 1.0   # Assumed to 1.0, not included in the calculations
dt = 0.01    # Timestep
nu = 1.0    # Kinematic viscosity
Lx = 1.0    # Canal length (x)
Ly = 1.0    # Canal length (y)
Nx = 10     # Quadratic elements (x)
Ny = 10     # Quadratic elements (y)

Vx0 = 0     # v(x,0).x
Vy0 = 0     # v(x,0).y
Vx1 = 0     # v(x,Ly).x
Vy1 = 0     # v(x,Ly).y

PA = 1.0    # Pressure at x=0
PB = 0.0    # Pressure at x=Lx

mesh = Rectangle(0,0,Lx,Ly,Nx,Ny,'left')

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V) # Trial function for u
p = TrialFunction(Q) # Trial function for p

v = TestFunction(V) # Test function for u
q = TestFunction(Q) # Test function for p

# Normal vector
n = FacetNormal(mesh)

p_in = Constant(PA) # Pressure at x=0
p_out = Constant(PB) # Pressure at x=1

def u0_boundary_bottom(x, on_boundary):
    if on_boundary:
        return x[1] < DOLFIN_EPS

def u0_boundary_top(x, on_boundary):
    if on_boundary:
        return x[1] > Ly - DOLFIN_EPS

def pressure_boundary_left(x):
    return x[0] < DOLFIN_EPS

def pressure_boundary_right(x):
    return x[0] > Lx - DOLFIN_EPS

# No slip makes sure that the velocity field is fixed at the boundaries
noslip_1  = DirichletBC(V, (0, 0), u0_boundary_bottom)
noslip_2  = DirichletBC(V, (0, 0), u0_boundary_top)

# Pressure boundary conditions
inflow  = DirichletBC(Q, p_in, pressure_boundary_left)
outflow = DirichletBC(Q, p_out, pressure_boundary_right)

bcu = [noslip_1, noslip_2]
bcp = [inflow, outflow]

# Create functions to save solutions
u0 = Function(V)
u1 = Function(V)
p0 = Function(Q)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

U = 0.5*(u0 + u)


# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx \
    + nu*inner(nabla_grad(U), nabla_grad(v))*dx \
    - p0*div(v)*dx\
    - inner(f, v)*dx\
    + inner(dot(u0,nabla_grad(u0)), v)*dx \
    + inner(p0*n,v)*ds\
    
a1 = lhs(F1)
L1 = rhs(F1)

# Calculation of phi
a2 = inner(grad(p), grad(q))*dx
L2 = inner(grad(p0),grad(q)) *dx\
    -(1.0/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1 - p0), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Time-stepping
t = 0.0

def step():
    global t
    # Compute tentative velocity step
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "amg")
    
    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", "amg")

    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "amg")

    # Move to next time step
    u0.assign(u1)
    p0.assign(p1)
    t += dt
    
    return t,u0,p1