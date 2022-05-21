#this script uses the bisection method to find a root of several continuous functions on a closed interval
#3.1.7C in the textbook


import numpy as np
import math
import matplotlib.pyplot as plt


def bisect(f,a,b,nmax,epsilon):
    a_0 = a
    b_0 = b
    n_thresh = (np.log(b-a)-np.log(2*epsilon))/np.log(2)
    fa = f(a)
    fb = f(b)
    if (fa > 0 and fb>0) or (fa<0 and fb<0):
        print('the function has the same signs at the endpoints of the interval, so the bisection method cannot be used to find any zeros contained in this interval')
    error = b - a 

    for i in range(0,nmax):
        error = error/2
        c = a + error
        fc = f(c)
        if fc == 0.0:
            print('we have found a zero at ',c)
            break
        print('iteration #',i+1)
        print('absolute error:','%.60f'%error)
        print('epsilon:       ','%.60f'%epsilon,'\n')
        
        if np.abs(error) < epsilon:
            print('The algorithm converged after', i+1,'iterations about the root x =',c,' of f(x) in the interval [',a_0,',',b_0,'] with an absolute error of at most',format(error,'.60g'))
            print('\nThe number of steps n required to reach an absolute error strictly less than this error tolerance using the bisection method is such that n > (log(b-a)-log(2*epsilon))/log(2)')
            print('so with parameters a =',a_0,', b =',b_0,', and epsilon =',epsilon,', n must be a positive integer strictly larger than',n_thresh,', which agrees with our computations.\n\n\n')
            break
        if fa<0 and fc>0:
            b=c
            fb=fc
        if fa>0 and fc<0:
            b=c
            fb=fc
        if fa<0 and fc<0:
            a=c
            fa=fc
        if fa>0 and fc>0:
            a=c
            fa=fc


def f(x):
    return pow(x,3)+3*x-1


def g(x):
    return pow(x,3)-2*np.sin(x)


def h(x):
    return x+10-x*np.cosh(50/x)


af = 0
bf = 1
ag = 0.5
bg = 2
ah = 120
bh = 130
nmax = 100
epsilon = 1.0/(2.0**24)

bisect(f,af,bf,nmax,epsilon)
bisect(g,ag,bg,nmax,epsilon)
bisect(h,ah,bh,nmax,epsilon)


x1 = np.linspace(0,1,100)
x2 = np.linspace(0.5,2,100)
x3 = np.linspace(120,130,100)
y1 = f(x1)
y2 = g(x2)
y3 = h(x3)

plt.plot(x1,y1)
plt.xlim(0,1)
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.title('f(x)')
plt.show()

plt.plot(x2,y2)
plt.xlim(0.5,2)
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.title('g(x)')
plt.show()

plt.plot(x3,y3)
plt.xlim(120,130)
plt.xlabel('x - axis')
plt.ylabel('y - axis')
plt.title('h(x)')
plt.show()

