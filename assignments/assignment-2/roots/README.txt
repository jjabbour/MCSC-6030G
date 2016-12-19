
The main program is in testroots.f90

The module newtonroots.f90 contains:

  function find_zero(x0) that implements Newton's method on
      the function f(x) = x - cos(x) in order to estimate the 
      zero of the function from a starting point x0.

  function fvals(xmin, xmax) that evaluates f(x) at
      nvals points between xmin and xmax and print them out to
      the file fvals.txt for  plotting purposes.
  
  parameters 
      maxiter: maximum number of iterations for Newton's method
      nvals:   number of points to evaluate f(x) for plotting.

The Python script makeplot.py creates a plot of f(x) using the
values in fvals.txt and the Newton iteration process using the 
values in newtonpts.txt, which are also written out as the iteration
proceeds.

To run the code:
  make test

To plot the results after running the code:
  make plot

The plot will be in plot.png

To clean up object code and executable:
  make clean

To also remove output files and the plot:
  make clobber

