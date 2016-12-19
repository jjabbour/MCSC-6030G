# initialize the function and its derivative intrested in for Newton solver
 
import numpy as np

def f(x):
    return x - np.cos(x)

def fprime(x):
    return 1 + np.sin(x)
