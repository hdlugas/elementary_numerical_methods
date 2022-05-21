#this script compares the secant method with Newton's method for a given differentiable function
#exercise 3.3.6C

import numpy as np
import matplotlib.pyplot as plt

def NewtonVsSecant(f,df,x,nmax,epsilon):
    fx = f(x)
    a = x
    for i in range(0,2):
        fp = df(x)
        d = fx/fp
        x = x-d
        if i == 1:
            x1 = x
        fx = f(x)

    b = x1
    fa = f(a)
    fb = f(b)
    if np.abs(fa) > np.abs(fb):
        tempx = a
        a = b
        b = tempx
        tempf = fa
        fa = fb
        fb = tempf


    for i in range(0,nmax-1):
        fp = df(x)
        d = fx/fp
        x = x-d
        fx = f(x)
            
        if i > 0:
            if np.abs(fa) > np.abs(fb):
                tempx = a
                a = b
                b = tempx
                tempf = fa
                fa = fb
                fb = tempf
            d_sec = (b-a)/(fb-fa)
            b = a
            fb = fa
            d_sec = d_sec*fa
            a = a-d_sec
            fa = f(a)
            print('\nSecant Method iteration: ',i)
            print('x =',x)
            print('f(x) =',fx)
            print('(b-a)/(f(b)-f(a)) =',d_sec)
            print('epsilon =',epsilon,'\n')
            if np.abs(d_sec) < epsilon:
                print('Since |(b-a)/(f(b)-f(a))| < epsilon, the Secant Method algorithm has converged at x =',a,'.\n')

        print('\nNewton Method iteration: ',i+1)
        print('x =',x)
        print('f(x) =',fx)
        print('f(x)/f\'(x) =',d)
        print('epsilon =',epsilon)
        if np.abs(d) < epsilon:
            print('Since f(x)/f\'(x) < epsilon, the Newton Method algorithm has converged at x =',x,'.\n')



def f(x):
    return pow(x,3)-3*x+1

def df(x):
    return 3*pow(x,2)-3

x0 = 2
nmax = 8
epsilon = 5*10**(-7)

NewtonVsSecant(f,df,x0,nmax,epsilon)


x = np.linspace(0,3,100)
f = np.zeros(100)
g = np.zeros(100)
for i in range(0,len(x)):
    f[i] = pow(x[i],3)-3*x[i]+1

plt.scatter(x,f,label='f(x)')
plt.scatter(x,g,label='g(x)=0')
plt.legend()
plt.show()
