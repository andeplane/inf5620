from ns_couette import *
from matplotlib import *


dt = 0.01 # Timestep
nu = 1.0 # Kinematic viscosity
rho = 1.0    # Density
#update_mesh(1,1,2,2)
#update_boundary_conditions(1,0,0,0)

t = 0.0
T = 4

#def vector_field_to_scalar_field(u):

V2 = VectorFunctionSpace(mesh, "CG", 1)

#while(t<T):
t,u = step()
coor = mesh.coordinates()
#print mesh.num_vertices()
#print len(u)
#print mesh

u_n = interpolate(u,V2)
u_vec = u_n.vector()
u_array = u_vec.array()

if 2*mesh.num_vertices() == len(u_array):
	for i in range(mesh.num_vertices()):
		print 'u(%8g,%8g) = %g' % (coor[i][0], coor[i,1], u_array[i])
