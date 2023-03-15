import numpy as np
def Calc_F_C_N(B : float,h : float,a_e : float,ecc2 : float) -> np.array:
    sinB = np.sin(B)
    sinB2 = sinB * sinB

    R_N = a_e * (1 - ecc2) / (1 - ecc2 * sinB2) ** 1.5 + h

    R_E = a_e / np.sqrt(1 - ecc2 * sinB2) + h

    F_C_N = np.array([
        [0.0,0.0,1.0 / R_E],
        [0.0,0.0,np.tan(B)/R_E],
        [-1.0/R_N,0.0,0.0]])

    return F_C_N