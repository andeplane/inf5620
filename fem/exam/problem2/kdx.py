from numpy import *
import matplotlib.pyplot as plt

case = [20.0,-1]
case = [0.5,-1]
case = [0.01,0.75]

C = case[0]
ymin = case[1]

alpha = 1.0
dx = 1.0

dt = C*dx**2/alpha
#C = alpha*dt/dx**2
print C

k = linspace(0,3.5,20)

a_exact = exp(-k**2*alpha*dt)

theta = 0
a_fe = (1 - 4*C*sin(k*dx/2)**2*(1-theta))/(1 + 4*C*sin(k*dx/2)**2*theta)
theta = 0.5
a_cn = (1 - 4*C*sin(k*dx/2)**2*(1-theta))/(1 + 4*C*sin(k*dx/2)**2*theta)
theta = 1
a_be = (1 - 4*C*sin(k*dx/2)**2*(1-theta))/(1 + 4*C*sin(k*dx/2)**2*theta)


plt.hold('on')
plt.plot(k*dx,a_exact,'bo-')
plt.plot(k*dx,a_fe,'ko-')
plt.plot(k*dx,a_be,'mo-')
plt.plot(k*dx,a_cn,'go-')
plt.legend(['Exact','Fe','BE','CN'])
title = "C=",C
plt.title(title)
plt.xlabel('k*dx')
plt.ylabel('Amplification factor')
plt.ylim([ymin,1])
#plt.xlim([0, 1])

plt.show()
