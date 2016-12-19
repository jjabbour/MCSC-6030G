# This algorithm solves the system of Poisson equations with homogeneous Dirichlet boundary conditions
# Delta u = v                   Delta v = f(x,y)
# v = 0 			v = 0 

# derived from the original problem the biharmonic equation on a Cartesian unit square   
# Delta^2 u = f(x,y)  on [0,1]x[0,1]
# u = 0 = Delta u
# The linear system of equation obtained from the discretisation is solved iteratively using Jacobi iteration     
# Based on oned.f from Using MPI (Gropp et. al.)

import sys
from mpi4py import MPI
import mpi_utils as mpiu
import numpy as np
import exchange
import matplotlib.pyplot as plt
from pylab import figure, show
from mpl_toolkits.mplot3d import Axes3D

comm = MPI.COMM_WORLD
NP = comm.Get_size()
pid = comm.Get_rank()

# Synchronize and start the timer
comm.Barrier()
local_time = np.array([-MPI.Wtime()],dtype = float)
start_time = np.array((0.0),dtype = float)

# Create a new 1D Cartesian communicator for a decomposition of the domain
comm1d = comm.Create_cart(dims=[NP], reorder=True)
# Get my position in this communicator, and my neighbors
myid = comm1d.Get_rank()
(left, right) = comm1d.Shift(direction=0, disp=1)

# Get problem size and compute parameters for block decomposition
# Initialise buffer to broadcast problem size
Nbuff = np.zeros(2, dtype=np.int)
Tbuff = np.zeros(1, dtype=np.float)
if pid == 0:
    n = 50
    TOL = 1.0e-6
    MAXITNS = 1000000
    nargin = len(sys.argv)
    if nargin>3:
        MAXITNS = int(sys.argv[3])
    if nargin>2:
        TOL = float(sys.argv[2])
    if nargin>1:
        n = int(sys.argv[1])
    Nbuff[0] = n
    Nbuff[1] = MAXITNS
    Tbuff[0] = TOL
comm.Bcast(Nbuff,root=0)
comm.Bcast(Tbuff,root=0)
nx = Nbuff[0]
ny = nx
MAXITNS = Nbuff[1]
TOL = Tbuff[0]
xs = np.linspace(0,1,nx+2)
ys = xs


lo_local = mpiu.BLOCK_LO(myid,nx,NP)+1
hi_local = mpiu.BLOCK_HI(myid,nx,NP)+1
nx_local = mpiu.BLOCK_SIZE(myid,nx,NP)

# Initialise the problem for computing solutions
# Problem setup: define RHS vector and mesh size
# f(x,y) is RHS when exact solution is u(x,y)= sin(m*pi*x)*sin(p*pi*y) on [0,1]x[0,1]
alpha = 1
beta  = 1
f = lambda x,y: (alpha**2+beta**2)**2*(np.pi)**4*np.sin(alpha*np.pi*x)*np.sin(beta*np.pi*y)

h = 1.0/(nx+1)
x = np.linspace(0,1,nx+2)
y = np.linspace(0,1,ny+2)
# Note reversal of arguments in call to meshgrid
(X,Y) = np.meshgrid(y,x[lo_local:hi_local+1]) # local grid for local grid points
# Notice: Y stores local (interior) y-coordinates only

# Observe: RHS vector must be scaled by h**2 for correct solution
# f must be sampled at interior points only.
F_RHS = (h**2)*f(X[:,1:-1],Y[:,1:-1])


# Now start actual iteration
def jacobi_sweep_2d(u,f):
    return 0.25*(f + u[:-2,1:-1] + u[2:,1:-1] + u[1:-1,:-2] + u[1:-1,2:]) - u[1:-1,1:-1]

def compute_norms(comm, delta_local, u_local):
    norm_delta_local = np.linalg.norm(delta_local.ravel())
    norm_u_local = np.linalg.norm(u_local.ravel())
    local_norms = np.array([norm_delta_local,norm_u_local])
    global_norms = np.zeros(2,dtype=np.float)
    comm.Allreduce([local_norms,2,MPI.DOUBLE],[global_norms,2,MPI.DOUBLE])
    global_norms = np.sqrt(global_norms)
    return global_norms

V = np.zeros((nx_local+2,ny+2))
U = np.zeros((nx_local+2,ny+2))
V[1:-1,1:-1] = 0.0
U[1:-1,1:-1] = 0.0
delta_U = np.zeros((nx_local,ny))
for itn in range(MAXITNS):
    exchange.exchange1d(comm1d,V,left,right)
    delta_V = jacobi_sweep_2d(V,F_RHS)
    V[1:-1,1:-1] +=  delta_V
    norms = compute_norms(comm,delta_V,V[1:-1,1:-1])
    converged1 = (norms[0]<=TOL*norms[1])
    
    exchange.exchange1d(comm1d,U,left,right)
    delta_U = jacobi_sweep_2d(U,V[1:-1,1:-1])
    U[1:-1,1:-1] +=  delta_U
    norms = compute_norms(comm,delta_U,U[1:-1,1:-1])
    converged2 = (norms[0]<=TOL*norms[1])
    converged = (converged1) and (converged2)
    if (converged):
        break

# Now assessing accuracy...
# Having obtained the solution, verify it
if not(converged):
    print "Jacobi iteration failed to converge after", itn, "iterations."
else:
    comm.Barrier()
    local_time[0] += MPI.Wtime()
    comm.Reduce(local_time,start_time,MPI.MAX,0)
    U_exact = np.sin(alpha*np.pi*X[:,1:-1])*np.sin(beta*np.pi*Y[:,1:-1])   
##### Write the results to an output file by one processor        
    if pid == 0:
        err = 1.0/n*(U_exact-U[1:-1,1:-1])
        outputfile = open('output.txt','a')
	outputfile.write(20*'-----'+'\n')
	outputfile.write('Results: \n')
	outputfile.write('-------- \n')	
	outputfile.write('The domain is discretised using '+str(n+2) + ' grid points. \n')
	outputfile.write('Jacobi iteration converges after '+ str(itn) + " iterations.\n")
        outputfile.write('The time of computation required is ' + str(start_time) + 'using '+ str(NP) + ' processes.' + '\n')

# cleaning up...
comm1d.Free() # Free 1D communicator
