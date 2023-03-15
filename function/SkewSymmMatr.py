import numpy as np
def SkewSymmMatr(a : np.array) -> np.array:

    Ax= np.array(
        [[0.0,-a[2],a[1]],
        [a[2],0.0,-a[0]],
        [-a[1],a[0],0.0]])
    return Ax