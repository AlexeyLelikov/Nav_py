import numpy
import numpy as np

x = numpy.ones((3,3))
print(x[:3,1].reshape((-1,1))+x[:3,1] *10)
print(x)
omega_IE_N_extrap = numpy.array([1, 2, 1]).reshape((-1,1))
print(omega_IE_N_extrap)

print(x[:3,1] + numpy.array([1, 2, 1]))


acc = np.array([1,2,3])
S = np.zeros((6,2))
S[0:3,1] = acc
print(S)

N = 10
queue_Lat = np.zeros(N)

print(queue_Lat)

