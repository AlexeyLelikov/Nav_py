import numpy as np
def DCM_nb(heading : float, pitch : float, roll : float) -> np.array:
    sz = np.sin(-heading)
    cz = np.cos(-heading)
    sy = np.sin(roll)
    cy = np.cos(roll)
    sx = np.sin(pitch)
    cx = np.cos(pitch)

    Cx = np.array([[1, 0, 0],[0, cx, sx],[0, -sx, cx]],dtype = np.float64)
    Cy = np.array([[cy, 0, -sy],[0, 1, 0],[sy, 0, cy]],dtype = np.float64)
    Cz = np.array([[cz, sz, 0],[-sz, cz, 0],[0, 0, 1]],dtype = np.float64)

    Cnb = Cy @ Cx @ Cz

    return Cnb

g = 9.81
acc_nav = np.array([0, 0, g]) #ENUp
alfa = 0
N = 6
C_N_B = np.eye(3)
acc_mas_x = np.array([])
acc_mas_y = np.array([])
for i in range(N+1):
    print(alfa * 180 / np.pi)
    acc_x = DCM_nb(0,alfa,0) @ acc_nav
    print("x",acc_x)
    print("norm",np.linalg.norm(acc_x))
    acc_y = DCM_nb(0,0,alfa) @ acc_nav
    alfa = alfa + 2 * np.pi / N
    print("y", acc_y)



