# $MYHG/assignments/assignment-2/roots/testzero.py
import newtonzero
print "Test routine for computing the zeros of function "
x0 = float(raw_input("Input x0, initial guess: "))
y = newtonzero.find_zero(x0)
print 'Estimated value: ', y

