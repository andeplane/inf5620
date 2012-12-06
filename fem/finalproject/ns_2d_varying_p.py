from dolfin import *
import numpy as np

#mesh = UnitSquare(10,10)
mesh = Rectangle(0,0,1,2,10,10,'left')

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

# Set parameter values
dt = 0.01 # Timestep
T = 40     # Final time
nu = 1.0 # Kinematic viscosity
rho = 1.0    # Density

#p_in = Constant(1.0) # Pressure at x=0
p_in = Expression("sin(2*pi*t)", t=0.0)
p_out = Constant(0.0) # Pressure at x=1
# No slip makes sure that the velocity field is zero at the boundaries
noslip_1  = DirichletBC(V, (0, 0),
                      "on_boundary && \
                       x[1] < DOLFIN_EPS")

noslip_2  = DirichletBC(V, (0, 0),
                      "on_boundary && \
                       x[1] > 2 - DOLFIN_EPS")

# We apply a pressure at x=0, no pressure at 
inflow  = DirichletBC(Q, p_in, "x[0] < DOLFIN_EPS")
outflow = DirichletBC(Q, p_out, "x[0] > 1.0 - DOLFIN_EPS")

bcu = [noslip_1, noslip_2]
bcp = [inflow, outflow]

# Create functions
u0 = Function(V)
u1 = Function(V)
p0 = Function(Q)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
U = 0.5*(u0 + u)

F1 = (1/k)*inner(u - u0, v)*dx \
    + nu*inner(grad(U), grad(v))*dx \
    - p0*div(v)*dx\
    - inner(f, v)*dx\
    + inner(grad(u0)*u0, v)*dx \
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

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:
    # Update pressure boundary condition
    p_in.t = t

    # Compute tentative velocity step
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    
    # Pressure correction
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", "amg")

    # Velocity correction
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")

    plot(u1, title="Velocity", rescale=True)

    # Save to file
    #ufile << u1
    #pfile << p1

    # Move to next time step
    u0.assign(u1)
    p0.assign(p1)
    t += dt

    u_vec = u1.vector()
    u_array = u_vec.array()
    print "t =", t
    print "max: ",u_array.max()


# Hold plot
interactive()
