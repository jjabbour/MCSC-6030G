import numpy as np 

def exp_taylor_np(x, N_terms):
    tol = 1.0e-16
    nmax = 100
    index = np.ones(len(x),dtype= int)
    term = np.ones(len(x),dtype= float)
    partial_sum = term.copy()
    mask = x < 0 
    x[mask] = -x 
    for j in range(1,N_terms+1):
	check = np.abs(term) > np.abs(partial_sum)*tol
    	index[check] = j
    	term *= np.where(check,x/float(j),1)
    	partial_sum += np.where(check,term,0)
    	if check.all == False: 
    	    break 
    partial_sum = np.where(mask,1/partial_sum,partial_sum)
    x[mask] = -x	
    return (partial_sum,index)       # returns value & number or terms summed

def test_exp_taylor_np():
    N_max = 100
    x = np.linspace(-20,20,11) 
    header = "x".center(10) + "# terms".rjust(8) + "true".rjust(19) 
    header += "approximate".rjust(27) + "error".rjust(21) 
    print header
    print "="*(10+12+24+24+24)
    (y,N_terms) = exp_taylor_np(x,N_max)
    exp_true = np.exp(x)
    relative_error = np.divide(abs(y-exp_true),exp_true)
    for k in range(0,len(x)-1):
        print repr(x[k]).center(10),repr(N_terms[k]).center(8),repr(exp_true[k]).rjust(24),repr(y[k]).rjust(24), \
	repr(relative_error[k]).rjust(24)




