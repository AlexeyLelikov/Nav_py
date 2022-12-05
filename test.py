import numpy

x = numpy.ones((3,3))
print(x[:3,1].reshape((-1,1))+x[:3,1] *10)
print(x)
omega_IE_N_extrap = numpy.array([1, 2, 1]).reshape((-1,1))
print(omega_IE_N_extrap)

print(x[:3,1] + numpy.array([1, 2, 1]))