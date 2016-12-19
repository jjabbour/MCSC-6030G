
import numpy as np

# jacobi_sweep_2d: solve the 2D poisson equation \Delta u = f(x,y) 
# with homogeneous Dirichlet boundary condition.
# the function returns the correction to the approximated solution. 

def jacobi_sweep_2d(u,f):
    return 0.25 * (f[1:-1,1:-1] + u[:-2,1:-1] + u[2:,1:-1]  \
                     + u[1:-1,:-2] + u[1:-1,2:]) - u[1:-1,1:-1]


# check_convergence: check the convergence of the jacobi iteration 
# based on the relative error of the 2 norm of the correction and the solution  
# |deltaU| / |U|  <  tolerance.

def check_convergence(diffU,U,tol):
    return np.linalg.norm(diffU.ravel()) <= np.linalg.norm(U.ravel())*tol


