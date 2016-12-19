# Solve the two dimensional Biharmonic equation \nabla^4 = f(x,y) on a unit square  
# with homogeneous Dirichelet boundary condition, \partial{u}\partial{n} = 0 on the boundary;

  
import sys
import numpy as np
import matplotlib.pyplot as plt
from pylab import figure, show 
from mpl_toolkits.mplot3d import Axes3D
import fun_biharmonic_jacobi as bj
# Default values for input parameters
n = 200
TOL = 1.0e-6
MAXITNS = 100000
nargin = len(sys.argv)
print "nargin =", nargin
if nargin>3:
   MAXITNS = int(sys.argv[3])
if nargin>2:
   TOL = float(sys.argv[2])
if nargin>1:
   n = int(sys.argv[1])

# Problem setup: define RHS vector and mesh size
a = 1.0
h = a/(n+1)
x = np.linspace(0,a,n+2)
(x,y) = np.meshgrid(x,x)
# f(x,y) is RHS when exact solution is u(x)= sin(m*pi*x)*sin(p*pi*y)
m = 1 
p = 1
f = lambda x,y: ((m/a)**2+(p/a)**2)**2*(np.pi**4)*np.sin(m*np.pi*x)*np.sin(p*np.pi*y)
# Observe: RHS vector must be scaled by h**2 for correct solution
#          f must be sampled at interior points only, 
F_RHS = np.zeros([n+2,n+2])
F_RHS[1:-1,1:-1] = (h**2)*f(x[1:-1,1:-1],y[1:-1,1:-1])
# for the given f(x,y), u(x,y) = sin(pi*x)*sin(pi*y) 
U_exact = np.sin(m*np.pi*x)*np.sin(p*np.pi*y)

# Solving the Biharmonic equation using Jacobi iteration 
U = np.zeros((n+2,n+2))
V = np.zeros((n+2,n+2))
delta_U = np.zeros([n,n])
delta_V = np.zeros([n,n])

for itn in range(MAXITNS):
    delta_V = bj.jacobi_sweep_2d(V,F_RHS)
    V[1:-1,1:-1] +=  delta_V
    delta_U = bj.jacobi_sweep_2d(U,h**2*V)
    U[1:-1,1:-1] +=  delta_U
    converged1 = bj.check_convergence(delta_U,U[1:-1,1:-1],TOL)
    converged2 = bj.check_convergence(delta_V,V[1:-1,1:-1],TOL)
    converged  = (converged1) and (converged2)
    if (converged):
        print itn 
        print converged
        break
   
# Having obtained the solution, verify it
if not(converged):
    print "Jacobi iteration failed to converge after", itn, "iterations."
else:
    err = (U_exact-U)
    print "Completed Jacobi iteration after", itn, "iterations."
    #print err.flatten()
    fig1 = figure()
    ax = Axes3D(fig1)
    plt.contourf(x,y,U,100)
    plt.colorbar()
    plt.show()    	


    print "The error between the exact and the numerical solution: ", \
    np.linalg.norm(err.flatten())
    fig2 = figure()
    plt.contourf(x,y,err,100)
    plt.colorbar()
    plt.show()

    fig3 = figure()
    ax = Axes3D(fig3)
    plt.contourf(x,y,U_exact,100)
    plt.colorbar()
    plt.show()

