#this script uses the central difference formula to approximate the derivative of a function

import numpy as np

def f(x):
    return np.sin(x)


x = np.pi/3
n = 6
h = 1
for i in range(1,n):
    dfdx = np.cos(x)
    dfdx1 = (f(x+h)-f(x-h))/(2*h)
    dfdx2 = (f(x+h/2)-f(x-h/2))/h
    error1 = np.abs(dfdx1-dfdx)
    error2 = np.abs(dfdx2-dfdx)
    print('\nIteration #:',i)
    print('h = ',h)
    print('\ndf/dx =',dfdx1,'at x =',x,'using h')
    print('df/dx =',dfdx2,'at x =',x,'using h/2')
    print('df/dx =',dfdx,'at x =',x,'is the actual value of df/dx at this point\n')
    print('absolute error using h:',error1)
    print('absolute error using h/2:',error2)
    print('\nThe ratio of the absolute error using h to the absolute error using h/2 is',error1/error2,'~ 4\n\n')
    h = h/2

#note that as expected, the error using h/2 is smaller than the error using h. As expected from the order of the error terms, the ratio of these terms is approxiately 4

