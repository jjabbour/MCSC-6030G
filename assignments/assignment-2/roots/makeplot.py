
import pylab
import numpy as np
import sys

# Read in the data files:
# Note: they might not exist, so use try...except:

try:
    x,fx = np.loadtxt('fvals.txt', unpack=True)
    xi,fi = np.loadtxt('newtonpts.txt', unpack=True)
except:
    # Reach this point if there's a problem reading the files
    print "Problem reading fvals.txt or newtonpts.txt, do they exist?"
    sys.exit(0)   # give up and quit!

pylab.figure(1)
pylab.clf()

# Plot the function f using the values from fvals.txt:
pylab.plot(x, fx, 'b')

# plot the Newton iteration process as red lines:
for i in range(len(xi)-1):
    pylab.plot([xi[i],xi[i],xi[i+1]], [0., fi[i], 0.], 'r')

# plot x-axis as a black line:
pylab.plot([x[0], x[-1]], [0., 0.], 'k')

# label x0, the initial guess:
pylab.plot([xi[0]], [0], 'k^')

# add an informative title:
pylab.title("Newton iteration process (triangle shows x0)")

# Save plot to a file:
fname = 'plot.png'
pylab.savefig(fname)
print "Created plot in file ",fname

