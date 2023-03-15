import numpy as np

def DCM_bn(heading : float, pitch : float, roll : float) -> np.array:
    sz = np.sin(pitch)
    cz = np.cos(pitch)
    sy = np.sin(-heading)
    cy = np.cos(-heading)
    sx = np.sin(roll)
    cx = np.cos(roll)

    Cx = np.array([[1, 0, 0],[0, cx, sx],[0, -sx, cx]],dtype = np.float64)
    Cy = np.array([[cy, 0, -sy],[0, 1, 0],[sy, 0, cy]],dtype = np.float64)
    Cz = np.array([[cz, sz, 0],[-sz, cz, 0],[0, 0, 1]],dtype = np.float64)

    Cnb = Cx @ Cz @ Cy

    return Cnb.T
