#solves a first-order differential equation
#with a given initial value using the Taylor Series Method of order 1,
#,i.e. Euler's Method
#exercise 7.1.C4


import numpy as np


def f(t,x):
    return x/(1+t)


def euler(a,b,n,x0,f,I):
    h = (b-a)/n
    t = a
    x = x0
    for i in range(0,n):
        x = x+h*f(t,x)
        t = t+h
        print('Iteration #:',i+1,'x(',t,') =',x)
        if i == n-1:
            print('Actual value of x(',1.0,')                =',I)
            print('\nAbsolute Error =',np.abs(I-x))


a = 0
b = 1#this is the value at which we desire the solution function to be evaluated at
n = 100
I = 2
x0 = 1

euler(a,b,n,x0,f,I)
