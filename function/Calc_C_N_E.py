import numpy as np
def Calc_C_N_E(B : float,L : float) -> np.array:

    cosL = np.cos(L)
    sinL = np.sin(L)
    cosB = np.cos(B)
    sinB = np.sin(B)

    C_N_E = np.zeros((3, 3))
    C_N_E[0, 0] = -cosL * sinB
    C_N_E[0, 1] = cosL * cosB
    C_N_E[0, 2] = -sinL
    C_N_E[1, 0] = -sinL * sinB
    C_N_E[1, 1] = sinL * cosB
    C_N_E[1, 2] = cosL
    C_N_E[2, 0] = cosB
    C_N_E[2, 1] = sinB
    C_N_E[2, 2] = 0

    return C_N_E
