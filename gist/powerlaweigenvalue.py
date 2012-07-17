from pylab import *

length = 5

#random matrix
mat = matrix(random((length,length)))
samplevec = ones(length)

powered = matrix_power(mat, 100)

print(dot(powered , samplevec))
