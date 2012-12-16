from ns_couette import *
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
ax.set_xlim(0, 2)
ax.grid()

T = 4

def run(data):
	t,u = step()
	u_array = interpolate(u,V2).vector().array()

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
		k = np.linspace(1,500,500)
		vx_e = Vx0*(1 - y0/Ly - np.sum(2.0/(k*np.pi)*np.sin(k*np.pi*y0/Ly)*np.exp(-k**2*np.pi**2*t*nu/Ly**2)))
		v_e[i,j,0] = vx_e
		v_e[i,j,1] = 0
		error = abs(u_array[n] - vx_e)
		if(error > 0.9):
			print "y0:",y0
			print "u_array[n]:",u_array[n]
			print "vx_e:",vx_e
		if error > max_error: max_error = error
		avg_error += error
	avg_error/=mesh.num_vertices()
	error = v_e - v

	print "t=%g, max error=%g" %(t,max_error)
	line.set_data(y, v[0,:,0])
	#line.set_data(y, error[0,:,0])
	return line,

ani = animation.FuncAnimation(fig, run, None, blit=False, interval=1,
    repeat=False)
#while True:
#	run(None)

plt.show()