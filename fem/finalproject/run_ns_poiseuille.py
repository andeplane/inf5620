from ns_poiseuille import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

V2 = VectorFunctionSpace(mesh, "CG", 1)

# Extract the mesh coordinates to calculate exact solution
Nx += 1 # There are Nx+1 vertices on Nx elements
Ny += 1 # There are Ny+1 vertices on Ny elements
coor = mesh.coordinates()
x = np.zeros(Nx)
y = np.zeros(Ny)
n=0
for l in range(Ny):
	y[l] = coor[n][1]
	for k in range(Nx):
		x[k] = coor[n][0]
		n += 1

#Prepare animation
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot([], [], lw=2)
ax.set_ylim(0, Vx0)
ax.set_xlim(0, 1)
ax.grid()

T = 4

def run(data):
	t,u = step()
	u_array = interpolate(u,V2).vector().array()
	u_array2 = interpolate(u,V).vector().array()

	avg_error = 0

	v_e = np.zeros( (Nx,Ny,2) )
	v = np.zeros( (Nx,Ny,2) )

	max_error = 0
	
	for n in range(mesh.num_vertices()):
		i = n % Nx
		j = n / Nx
		x0 = coor[n][0]
		y0 = coor[n][1]
		v[i,j,0] = u_array[n]
		v[i,j,1] = u_array[n+mesh.num_vertices()]
		k = np.linspace(1,50,50)
		vx_e = (PA-PB)/(2*nu*rho*Lx)*y0*(Ly-y0)

		v_e[i,j,0] = vx_e
		v_e[i,j,1] = 0
		error = abs(u_array[n] - vx_e)
		if error > max_error: max_error = error
		avg_error += error
	avg_error/=mesh.num_vertices()
	error = np.log(abs(v_e - v))
	v_max = np.max(v)
	print "v_max: %g" % (v_max)

	print "t=%g, max error=%g" %(t,max_error)
	line.set_data(y, v[0,:,0])
	return line,

ani = animation.FuncAnimation(fig, run, None, blit=False, interval=1,
    repeat=False)
#while True:
#	run(None)

plt.show()