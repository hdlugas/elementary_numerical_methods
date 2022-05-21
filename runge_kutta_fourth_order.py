#this script uses the fourth-order Runge-Kutta method
#to approximate a first-order differential equation


import numpy as np
import matplotlib.pyplot as plt


def f(t,x):
    return 1/(x**2)-x*t


a = 1
b = 2
n = 2**6
I = 2
x0 = 1

def rk4(a,b,n,x0,f):
    h = (b-a)/n
    t = a
    x = x0
    t_vec = []
    x_vec = []
    for i in range(0,n):
        k1 = h*f(t,x)
        k2 = h*f(t+0.5*h,x+0.5*k1)
        k3 = h*f(t+0.5*h,x+0.5*k2)
        k4 = h*f(t+h,x+k3)
        x = x+(k1+2*k2+2*k3+k4)/6
        t = a+h*(i+1)
        print('\nstep #:',i+1)
        print('x =',x,'at t =',t)
        t_vec.append(t)
        x_vec.append(x)
    return(t_vec,x_vec)


t,x = rk4(a,b,n,x0,f)
plt.plot(t,x)
plt.xlim(1,2)
plt.xlabel('t -axis')
plt.ylabel('x-axis')
plt.title('Solution function f(x)')
plt.show()



