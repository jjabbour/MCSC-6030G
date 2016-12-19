# $CLASSHG/codes/python/roots/newtonroots.py
import functions as fun
maxiter = 50  # max number of Newton iterations
nvals = 1001  # number of points for plotting f
TOL = 1.e-15  # convergence tolerance

def find_zero(x0):
    ''' Estimate the n'th root of a using Newton's method. 
    Requires an initial guess x0
    Returns:
        the estimate, if Newton's method converged,
        999999 if it did not converge in maxiter iterations.
    '''
    x = x0 # initial guess
    # Keep track of the range of x values used in iteration,
    # for plotting purposes.  Initialize to initial point:
    xmin = x
    xmax = x
    # open file for writing iterates
    outfile = open('newtonpts.txt','w')

    # Newton iteration to find a zero of f(x) = x-cos(x).

    converged = False
    for j in range(2*maxiter+1):
        # evaluate function and its derivative:
        fx = fun.f(x)
        fxprime = fun.fprime(x)

        # compute Newton increment x:
        delta = fx/fxprime

        # print values to watch convergence:
        print "j=%3d, x=%25.15g, f(x)=%16.6g" % (j,x,fx)
        # save x and fx for plotting purposes:
        outfile.write('%16.6e %16.6e\n' % (x,fx))
        # update x:
        x -= delta 
        
        # update min and max seen so far:
        xmin = min(x, xmin)
        xmax = max(x, xmax)

        # check for convergence:
        if (abs(delta) < TOL):
            converged = True
            break

    outfile.close()
    print 'Number of iterations: %d' % j

    # Increase the range of x a bit for better plots:
    xtra = 0.2*(xmax - xmin)
    xmin -= xtra
    xmax += xtra

    # Print out values of x over expanded interval for plotting:
    fvals(xmin,xmax)
    # If iteration didn't converge, set a special value:
    if converged:
        return x
    else:
        return 999999.e0


#-------------------------------------------------------------------

def fvals(xmin, xmax):
    ''' Print out values of the function f(x) = x**n - a for plotting purposes,
    in the interval xmin to xmax.  
    nval, the number of values to print, is a module parameter.
    '''
    import numpy as np
    x = np.linspace(xmin,xmax,nvals)
    fx = fun.f(x)
    A = np.vstack([x, fx]).T
    np.savetxt('fvals.txt', A)
    return None
