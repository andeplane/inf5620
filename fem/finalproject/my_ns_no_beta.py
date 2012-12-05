from dolfin import *

mesh = UnitSquare(10,10)

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# Define trial and test functions
u = TrialFunction(V) # Trial function for u
p = TrialFunction(Q) # Trial function for p
phi = TrialFunction(Q) # Trial function for phi

v = TestFunction(V) # Test function for u
q = TestFunction(Q) # Test function for p

# Set parameter values
dt = 0.001 # Timestep
T = 4     # Final time
nu = 0.5 # Kinematic viscosity
rho = 1.0    # Density

p_in = Constant(1.0) # Pressure at x=0
#p_in = Expression("sin(2*pi*t)", t=0.0)
p_out = Constant(0.0) # Pressure at x=1
# No slip makes sure that the velocity field is zero at the boundaries
noslip_1  = DirichletBC(V, (0, 0),
                      "on_boundary && \
                       x[1] < DOLFIN_EPS")

noslip_2  = DirichletBC(V, (0, 0),
                      "on_boundary && \
                       x[1] > 1 - DOLFIN_EPS")

# We apply a pressure at x=0, no pressure at 
inflow  = DirichletBC(Q, p_in, "x[0] < DOLFIN_EPS")
outflow = DirichletBC(Q, p_out, "x[0] > 1.0 - DOLFIN_EPS")

bcu = [noslip_1, noslip_2]
bcp = [inflow, outflow]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx

a1 = lhs(F1)
L1 = rhs(F1)

# Calculation of pressure
a2 = inner(grad(p), grad(q))*dx
L2 = -(rho/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k/rho*inner(grad(p1), v)*dx

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
    #begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    #end()

    # Pressure correction
    #begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", "amg")
    #end()

    # Velocity correction
    #begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    #end()

    plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt
    print "t =", t

# Hold plot
interactive()
