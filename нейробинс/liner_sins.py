import torch
import math as m
import matplotlib.pyplot as plt

def DCM_bn(heading : float, pitch : float, roll : float) -> torch.Tensor:
    sz = m.sin(pitch)
    cz = m.cos(pitch)
    sy = m.sin(-heading)
    cy = m.cos(-heading)
    sx = m.sin(roll)
    cx = m.cos(roll)

    Cx = torch.Tensor([[1, 0, 0],[0, cx, sx],[0, -sx, cx]])
    Cy = torch.Tensor([[cy, 0, -sy],[0, 1, 0],[sy, 0, cy]])
    Cz = torch.Tensor([[cz, sz, 0],[-sz, cz, 0],[0, 0, 1]])

    Cnb = Cx @ Cz @ Cy

    return Cnb.T

def threeaxisrot(r11, r12, r21, r31, r32):
    r1 = m.atan2(r11,r12)
    r2 = m.asin(r21)
    r3 = m.atan2(r31,r32)
    return r1,r2,r3

def EulerFromDCMnue(dcm):
    (r1, r2, r3) = threeaxisrot(dcm[0, 1],dcm[0,0], -dcm[0, 2],dcm[1, 2], dcm[2, 2])
    return (r1, r2, r3)

def Skew(a : torch.Tensor) -> torch.Tensor:

    Ax= torch.Tensor(
        [[0.0,-a[2],a[1]],
        [a[2],0.0,-a[0]],
        [-a[1],a[0],0.0]])
    return Ax

def Calc_F_C_N(B : float, h : float, a_e : float, ecc2 : float) -> torch.Tensor:
    sinB = m.sin(B)
    sinB2 = sinB * sinB

    R_N = a_e * (1 - ecc2) / (1 - ecc2 * sinB2) ** 1.5 + h

    R_E = a_e / m.sqrt(1 - ecc2 * sinB2) + h

    F_C_N = torch.Tensor([
        [0.0,0.0,1.0 / R_E],
        [0.0,0.0,m.tan(B)/R_E],
        [-1.0/R_N,0.0,0.0]])

    return F_C_N
# огроничение для курса
def twoPiBound(heading : float) -> float:
    if heading < 0:
        heading = heading + 2.0 * m.pi
    return heading

def SINS(NavState : torch.Tensor, Sensors : torch.Tensor, C_i_b :  torch.Tensor, C_i_n : torch.Tensor, dt : float):
    # константы
    Om_e = 7.292115e-5
    a_e = 6378136.0
    ecc = 0.0818191065283638
    ecc2 = ecc * ecc
    R = 6356863 # м
    # параметры
    delta_v = Sensors[0:3]
    alfa = Sensors[3:6]
    B = NavState[0]
    L = NavState[1]
    Alt = NavState[2]
    W = NavState[3:6]
    sinB = m.sin(B)
    cosB = m.cos(B)
    sinL = m.sin(L)
    cosL = m.cos(L)
    sinB2 = sinB * sinB
    E = torch.Tensor([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])

    omega_IE_N = torch.Tensor([Om_e * cosB, Om_e * sinB, 0.0])
    omega_sum = omega_IE_N + torch.Tensor([ W[2] / R, W[2] / R * m.tan(B), -W[0] / R])
    skew_omega = Skew(omega_sum)

    skew_alfa = Skew(alfa)  # Коссосиметричная матрица
    mod_alpha = m.sqrt(alfa[0] * alfa[0] + alfa[1] * alfa[1] + alfa[2] * alfa[2])
    mod_alpha2 = mod_alpha * mod_alpha

    C_b_b = E - (1 - mod_alpha2 / 6) * skew_alfa + (0.5 - mod_alpha2 / 24) * skew_alfa @ skew_alfa
    C_i_b = C_b_b @ C_i_b  # Получили МНК из body в инерциальную в данный такт
    #C_B_N_d1 = E + Skew(gyro)
    C_n_n = E - dt * skew_omega + ((dt ** 2) / 2) * skew_omega @ skew_omega

    C_i_n = C_n_n @ C_i_n
    C_b_n = C_i_n @ C_i_b.T

    delta_W_n = C_b_n @ delta_v
    omega = omega_sum
    W_N = W[0] + delta_W_n[0] + dt * (-(omega_IE_N[1] + omega[1]) * W[2] + omega[2] * W[1])
    W_U = W[1] + delta_W_n[1] + dt * ((omega_IE_N[0] + omega[0]) * W[2] - omega[2] * W[0] - gt)
    W_E = W[2] + delta_W_n[2] + dt * ((omega_IE_N[1] + omega[1]) * W[0] - (omega_IE_N[0] + omega[0]) * W[1])

    Lat = B + dt * W_N / R  # Широта
    Lon = L + dt * W_E / (R * cosB)  # долгота

    roll = m.atan2(-C_b_n[1, 2], C_b_n[1, 1])
    pitch = m.atan(C_b_n[1, 0] / m.sqrt(C_b_n[1, 1] * C_b_n[1, 1] + C_b_n[1, 2] * C_b_n[1, 2]))
    heading = m.atan2(C_b_n[2, 0], C_b_n[0, 0])
    heading = twoPiBound(heading)

    NavState = torch.Tensor([Lat, Lon, Alt,  # Широта, долгота, геодезическая высота
                             W_N, W_U, W_E,  # Скорость на восток, север, вверх
                             roll, pitch, heading])  # Крен , тангаж , курс
    
    return NavState , C_i_b , C_i_n

Lat = 55.0 * (m.pi / 180.0)
Lon = 38.0 * (m.pi / 180.0)
Alt = 0

gt = 9.81
U = 7.292115e-5

bias_gyro_x = (0.01 * m.pi / 180) / 3600.0 # рад / c
bias_gyro_y = (0.01 * m.pi / 180) / 3600.0 # рад / c
bias_gyro_z = (0.01 * m.pi / 180) / 3600.0 # рад / c
bias_acc_x = 0.2 * 1e-3 * gt
bias_acc_y = 0.2 * 1e-3 * gt
bias_acc_z = 0.2 * 1e-3 * gt
# вектора смещений
bias_gyro = torch.Tensor([bias_gyro_y, bias_gyro_z, bias_gyro_x])
bias_acc = torch.Tensor([bias_acc_y, bias_acc_z, bias_acc_x])
E = torch.Tensor([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
K_gyro = torch.eye(3)
K_acc = torch.eye(3)
# идеальные показания
acc = torch.Tensor([0.0, gt, 0.0])
gyro = torch.Tensor([(U * m.cos(Lat)), (U * m.sin(Lat)), 0])
# Синтезированные показания
acc = K_acc @ acc + bias_acc
gyro = K_gyro @ gyro + bias_gyro
Sensors = torch.zeros((6,))
W_NUE = torch.Tensor([0,0,0])
dt = 1.0 / 100.0
Roll = 1 * m.pi / 180
Pitch = 2 * m.pi / 180
Heading = 180.0 * m.pi / 180
C_B_N = DCM_bn(Heading, Pitch, Roll)
roll = m.atan2(-C_B_N[1, 2], C_B_N[1, 1])
pitch = m.atan(C_B_N[1, 0] / m.sqrt(C_B_N[1, 1] * C_B_N[1, 1] + C_B_N[1, 2] * C_B_N[1, 2]))
heading = m.atan2(C_B_N[2, 0], C_B_N[0, 0])
heading = twoPiBound(heading)
print(heading * 180 / m.pi, pitch * 180 / m.pi, roll * 180 / m.pi)
C_I_B = C_B_N.T
C_I_N = E
acc = C_B_N.T @ acc
gyro = C_B_N.T @ gyro
T = 60 * 90 # c
N = 100 * T
t = torch.linspace(0,T,N)

NavState = torch.Tensor([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_NUE[0],W_NUE[1],W_NUE[2],
                    Roll,Pitch,Heading])           # Скорость на север, вверх, на восток

queue_Lat = torch.zeros(N)
queue_Lon = torch.zeros(N)
queue_pitch = torch.zeros(N)
queue_roll = torch.zeros(N)
queue_heading = torch.zeros(N)
queue_Alt = torch.zeros(N)
queue_W_N = torch.zeros(N)
queue_W_U = torch.zeros(N)
queue_W_E = torch.zeros(N)

Sensors[0:3] = acc * 1e-2 # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
Sensors[3:6] = gyro * 1e-2

for i in range(0,N):
    NavState[2] = 0.0
    NavState[4] = 0.0
    queue_Lat[i] = NavState[0]
    queue_Lon[i] = NavState[1]
    queue_roll[i] = NavState[6]
    queue_pitch[i] = NavState[7]
    queue_heading[i] = NavState[8]
    queue_W_N[i] = NavState[3]
    queue_W_U[i] = NavState[4]
    queue_W_E[i] = NavState[5]
    queue_Alt[i] = NavState[2]
    (NavState , C_I_B , C_I_N) = SINS(NavState, Sensors, C_I_B , C_I_N, dt)

queue_Lat = queue_Lat * (180.0 / m.pi)
queue_Lon = queue_Lon * (180.0 / m.pi)
queue_pitch = queue_pitch * (180.0 / m.pi)
queue_roll = queue_roll * (180.0 / m.pi)
queue_heading = queue_heading * (180.0 / m.pi)

plt.figure(1)
plt.plot(t,queue_Lat)
plt.title("Широта")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(2)
plt.plot(t,queue_Lon)
plt.title("Долгота")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(3)
plt.plot(t,queue_pitch)
plt.title("Тангаж")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(4)
plt.plot(t,queue_roll)
plt.title("Крен")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(5)
plt.plot(t,queue_heading)
plt.title("Курс")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(6)
plt.plot(t,queue_W_N)
plt.title("Скорость W_N")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

plt.figure(7)
plt.plot(t,queue_W_U)
plt.title("Скорость U")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

plt.figure(8)
plt.plot(t,queue_W_E)
plt.title("Скорость E")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

R = 6356.863 # км
queue_x = (queue_Lat * ((m.pi / 180) * R)) - (queue_Lat[0] * ((m.pi / 180) * R))
queue_y = (queue_Lon * (m.pi / 180) * (R * m.cos(Lat))) - (queue_Lon[0] * m.pi / 180 * (R * m.cos(Lat)))
plt.figure(9)
plt.plot(t,queue_x)
plt.title("N")
plt.xlabel('Cекунды')
plt.ylabel('N,км')
plt.grid(True)
plt.show()

plt.figure(10)
plt.plot(t,queue_y)
plt.title("E")
plt.xlabel('Cекунды')
plt.ylabel('E,км')
plt.grid(True)
plt.show()

print("done")