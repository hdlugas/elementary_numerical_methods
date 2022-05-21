#this script uses a finite difference method to approximate the solution of a parabolic PDE

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def parabolic_pde_solver_explicit(h,k,n,m):
    u = np.zeros((n+1))
    v = np.zeros((n+1))
    w = np.zeros((n+1))
    sigma = k/pow(h,2)
    sol_approx = np.zeros((m+1,n+1))#top row is solution at t_0, then time step increases as you move down columns
    sol_exact = np.zeros((m+1,n+1))
    for i in range(1,n):
        u[i] = np.sin(np.pi*i*h)
    sol_approx[0,:] = u
    sol_exact[0,:] = u
    for j in range(1,m+1):
        for i in range(1,n):
            v[i] = sigma*u[i+1]+(1-2*sigma)*u[i]+sigma*u[i-1]
            #v[i] = 0.5*(u[i-1]+u[i+1])
        sol_approx[j,:] = v
        t = j*k
        for i in range(1,n):
            w[i] = pow(np.e,-1*np.pi*np.pi*t)*np.sin(np.pi*i*h)
        sol_exact[j,:] = w
        for i in range(1,n):
            u[i] = v[i]
    return sol_approx,sol_exact


h = pow(2,-4)
k = pow(2,-10)

n = int((1/h))

m = int((0.25/k))


sol_approx,sol_exact = parabolic_pde_solver_explicit(h,k,n,m)
x = np.linspace(0,1,n+1)
y = np.linspace(0,0.25,m+1)
x,y = np.meshgrid(x,y)


fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_approx = ax.plot_surface(x,y,sol_approx,linewidth=0,antialiased=False,color='blue')
ax.set_title('Approximate Solution Surface')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()


fig,ax = plt.subplots(subplot_kw={"projection": "3d"})
surf_exact = ax.plot_surface(x,y,sol_exact,linewidth=0,antialiased=False,color='blue')
ax.set_title('Exact Solution Surface')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()


fig,ax = plt.subplots(1,1)
ax.contour(x,y,sol_approx,levels=50)
ax.contour(x,y,sol_exact,levels=50)
ax.set_title('Contour Plot of Exact vs Approximate Solution')
ax.set_xlabel('Space')
ax.set_ylabel('Time')
plt.show()
