import numpy as np
import matplotlib.pyplot as plt
from Savage import Savage
from LatLonAlt2XYZ import LatLonAlt2XYZ
from CalcGravity import CalcGravity
from DCM_bn import DCM_bn
from cross import cross
from SkewSymmMatr import SkewSymmMatr
from Calc_F_C_N import Calc_F_C_N
from twoPiBound import twoPiBound

Lat = 56.0 * (np.pi / 180.0)
Lon = 38.0 * (np.pi / 180.0)
Alt = 0


# Вычисление силы тяжести

(X,Y,Z) = LatLonAlt2XYZ(Lat, Lon, Alt)
gravity = CalcGravity(X, Y, Z)
sinB = np.sin(Lat)
cosB = np.cos(Lat)
sinL = np.sin(Lon)
cosL = np.cos(Lon)

g_P_N_extrap = np.zeros(3)
g_P_N_extrap[0] = -cosL * sinB * gravity[0] - sinL * sinB * gravity[1] + cosB * gravity[2]
g_P_N_extrap[1] = cosL * cosB * gravity[0] + sinL * cosB * gravity[1] + sinB * gravity[2]
gt = abs(g_P_N_extrap[1])
U = 7.292115e-5

bias_gyro_x = (30 * np.pi / 180) / 3600.0 # рад / c
bias_gyro_y = (30 * np.pi / 180) / 3600.0 # рад / c
bias_gyro_z = (30 * np.pi / 180) / 3600.0 # рад / c
bias_acc_x = 20 * 1e-3 * gt
bias_acc_y = 20 * 1e-3 * gt
bias_acc_z = 20 * 1e-3 * gt

# вектора смещений
bias_gyro = np.array([bias_gyro_y, bias_gyro_z, bias_gyro_x],dtype = np.float64)
bias_acc = np.array([bias_acc_y, bias_acc_z, bias_acc_x],dtype = np.float64)
K_gyro = np.eye(3,dtype = np.float64)
K_acc = np.eye(3,dtype = np.float64)
# идеальные показания
acc = np.array([0.0 , gt , 0.0 ], dtype = np.float64)
gyro = np.array([(U * np.cos(Lat)), (U * np.sin(Lat)) , 0],dtype = np.float64)
# Синтезированные показания
acc = K_acc @ acc + bias_acc
gyro = K_gyro @ gyro + bias_gyro

Sensors = np.zeros((6,),dtype = np.float64)
W_NUE = np.array([0,0,0])
W_NUE_old = np.array([0,0,0])
dt = 1.0 / 100.0
Roll = 2 * np.pi / 180
Pitch = 3 * np.pi / 180
Heading = 180.0 * np.pi / 180
C_B_N = DCM_bn(Heading, Pitch, Roll)
acc = C_B_N.T @ acc
gyro = C_B_N.T @ gyro
T = 60 * 90 # c
N = 100 * T
t = np.linspace(0,T,N)

NavState = np.array([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_NUE[0],W_NUE[1],W_NUE[2],
                    Roll,Pitch,Heading])

queue_Lat = np.zeros(N)
queue_Lon = np.zeros(N)
queue_pitch = np.zeros(N)
queue_roll = np.zeros(N)
queue_yaw = np.zeros(N)
queue_Alt = np.zeros(N)
queue_W_N = np.zeros(N)
queue_W_U = np.zeros(N)
queue_W_E = np.zeros(N)

Sensors[0:3] = acc * 1e-2 # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
Sensors[3:6] = gyro *1e-2
h4 = 1 / 25

# Статические переменные
Om_e = 7.292115e-5
E = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]],dtype = np.float64)
# Используется ПЗ-90
a_e = 6378136.0
ecc = 0.0818191065283638
ecc2 = ecc * ecc
iter = 0
delta_v_array = np.zeros((3,4))
delta_alfa_array = np.zeros((3,4))

C_I_B = C_B_N.T
C_I_N = E

def SINS(NavState : np.array, Sensors : np.array, dt : float, C_i_b : np.array, C_i_n : np.array , gt : float) -> tuple:

    global Om_e
    global E
    # Используется ПЗ-90
    global a_e
    global ecc2
    global iter
    global delta_v_array
    global delta_alfa_array

    # переменные квази-скоростей и квази-угллов
    delta_v = Sensors[0:3]  # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
    delta_alpha = Sensors[3:6]  # выходной сигнал гироскопов (скорость) по оси (x,y,z) на текущем такте расчетов
    delta_v_array[0:3, iter] = delta_v
    delta_alfa_array[0:3, iter] = delta_alpha

    if iter == 3:

        B = NavState[0]
        L = NavState[1]  # текущая  долгота ВС в выбранной СК
        h = NavState[2]  # текущая  геодезическая высота ВС относительно эллипсоида выбранной СК
        W = NavState[3:6]  # текущий вектор земной скорости ВС (N,U,E)
        sinB = np.sin(B)
        cosB = np.cos(B)
        sinL = np.sin(L)
        cosL = np.cos(L)
        sinB2 = sinB * sinB
        R = 6356863 # м
        # Расчет вектора конечного поворота Эйлера (body)
        alfa_1_2 = np.array([delta_alfa_array[0, 0], delta_alfa_array[1, 0], delta_alfa_array[2, 0]]) + np.array(
            [delta_alfa_array[0, 1], delta_alfa_array[1, 1], delta_alfa_array[2, 1]])
        alfa_3_4 = np.array([delta_alfa_array[0, 2], delta_alfa_array[1, 2], delta_alfa_array[2, 2]]) + np.array(
            [delta_alfa_array[0, 3], delta_alfa_array[1, 3], delta_alfa_array[2, 3]])
        alfa = alfa_1_2 + alfa_3_4
        beta = 2 / 3 * cross(alfa_1_2, alfa_3_4)  # коннинг
        alfa = alfa + beta  # вектор порота body эйлера
        skew_alfa = SkewSymmMatr(alfa)  # Коссосиметричная матрица
        mod_alpha = np.sqrt(alfa[0] * alfa[0] + alfa[1] * alfa[1] + alfa[2] * alfa[2])
        mod_alpha2 = mod_alpha * mod_alpha
        # Расчет вектора конечного поворота (скорости) Эйлера (Navigation)
        omega_IE_N = np.array([Om_e * cosB, Om_e * sinB, 0.0])
        F_C_N = Calc_F_C_N(B, h, a_e, ecc2)
        omega_sum = omega_IE_N + np.array([ W[2] / R, W[2] / R * np.tan(B), -W[0]/ R])
        skew_omega = SkewSymmMatr(omega_sum)
        # нахождение МНК из body t-1 в body t при помощи решения уравенения Пуасона
        C_b_b = E - (1 - mod_alpha2 / 6) * skew_alfa + (0.5 - mod_alpha2 / 24) * skew_alfa @ skew_alfa
        C_b_b = (E - ((1.0 / 2.0) * ( C_b_b @  C_b_b.T - E))) @  C_b_b
        C_i_b = C_b_b @ C_i_b  # Получили МНК из body в инерциальную в данный такт
        # нахождение МНК из навигационной СК t-1 в навигационную СК t при помощи решения уравенения Пуасона
        C_n_n = E - dt * skew_omega + ((dt ** 2) / 2) * skew_omega @ skew_omega
        C_n_n = ( E - ((1.0 / 2.0) * ( C_n_n @  C_n_n.T - E))) @  C_n_n
        C_i_n = C_n_n @ C_i_n
        C_b_n = C_i_n @ C_i_b.T
        # нормализация и ортогонализация мнк
        C_i_b = (E - ((1.0 / 2.0) * ( C_i_b @  C_i_b.T - E))) @  C_i_b
        C_i_n = (E - ((1.0 / 2.0) * ( C_i_n @  C_i_n.T - E))) @  C_i_n
        C_b_n = (E - ((1.0 / 2.0) * ( C_b_n @  C_b_n.T - E))) @  C_b_n
        # нахождение приращения линейных скоростей методом Рунге-Кутта 4-го порядка
        delta_W = np.array([0, 0, 0])
        h1 = 1 / 100
        for i in range(4):
            skew_delta_alfa = SkewSymmMatr(delta_alfa_array[0:3, i])
            w_b = delta_v_array[0:3, i]
            k_1 = w_b - skew_delta_alfa @ delta_W
            k_2 = w_b - skew_delta_alfa @ (delta_W + (h1 / 2) * k_1)
            k_3 = w_b - skew_delta_alfa @ (delta_W + (h1 / 2) * k_2)
            k_4 = w_b - skew_delta_alfa @ (delta_W + h1 * k_3)
            delta_W = delta_W + (k_1 + 2 * k_2 + 2 * k_3 + k_4) / 6

        delta_W_n = C_b_n @ delta_W

        omega = omega_sum
        W_N = W[0] + delta_W_n[0] + dt * (-(omega_IE_N[1] + omega[1]) * W[2] + omega[2] * W[1])
        W_U = W[1] + delta_W_n[1] + dt * ((omega_IE_N[0] + omega[0]) * W[2] - omega[2] * W[0] - gt)
        W_E = W[2] + delta_W_n[2] + dt * ((omega_IE_N[1] + omega[1]) * W[0] - (omega_IE_N[0] + omega[0]) * W[1])

        Lat = B + dt * W_N / R  # Широта
        Lon = L + dt * W_E / (R * cosB)  # долгота
        Alt = 0
        roll = np.arctan2(-C_b_n[1, 2], C_b_n[1, 1])
        pitch = np.arctan(C_b_n[1, 0] / np.sqrt(C_b_n[1, 1] * C_b_n[1, 1] + C_b_n[1, 2] * C_b_n[1, 2]))
        heading = np.arctan2(C_b_n[2, 0], C_b_n[0, 0])
        heading = twoPiBound(heading)  # Ограничение

        NavState = np.array([Lat, Lon, Alt,  # Широта, долгота, геодезическая высота
                             W_N, W_U, W_E,  # Скорость на север, вверх, на восток
                             roll, pitch, heading])  # Крен , тангаж , курс

        delta_v_array = np.zeros((3, 4))
        delta_alfa_array = np.zeros((3, 4))
        iter = 0
    else:
        iter += 1
    return NavState, C_i_b, C_i_n

for i in range(0,N):
    NavState[2] = 0.0
    NavState[4] = 0.0
    queue_Lat[i] = NavState[0]
    queue_Lon[i] = NavState[1]
    queue_roll[i] = NavState[6]
    queue_pitch[i] = NavState[7]
    queue_yaw[i] = NavState[8]
    queue_W_N[i] = NavState[3]
    queue_W_U[i] = NavState[4]
    queue_W_E[i] = NavState[5]
    queue_Alt[i] = NavState[2]
    (NavState, C_I_B, C_I_N) = SINS(NavState, Sensors, h4, C_I_B, C_I_N, gt)

queue_Lat = queue_Lat * (180.0 / np.pi)
queue_Lon = queue_Lon * (180.0 / np.pi)
queue_pitch = queue_pitch * (180.0 / np.pi)
queue_roll = queue_roll * (180.0 / np.pi)
queue_yaw = queue_yaw * (180.0 / np.pi)

plt.figure(1)
plt.plot(t,queue_Lat)
plt.title("Широта")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/1.png')
plt.show()

plt.figure(2)
plt.plot(t,queue_Lon)
plt.title("Долгота")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/2.png')
plt.show()

plt.figure(3)
plt.plot(t,queue_pitch)
plt.title("Тангаж")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/3.png')
plt.show()

plt.figure(4)
plt.plot(t,queue_roll)
plt.title("Крен")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/4.png')
plt.show()

plt.figure(5)
plt.plot(t,queue_yaw)
plt.title("Курс")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/5.png')
plt.show()

plt.figure(6)
plt.plot(t,queue_W_N)
plt.title("Скорость W")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.savefig('img/6.png')
plt.show()

plt.figure(7)
plt.plot(t,queue_W_U)
plt.title("Скорость N")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.savefig('img/7.png')
plt.show()

plt.figure(8)
plt.plot(t,queue_W_E)
plt.title("Скорость E")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.savefig('img/8.png')
plt.show()

R = 6356.863 # км
queue_x = (queue_Lat * ((np.pi / 180) * R)) - (queue_Lat[0] * ((np.pi / 180) * R))
queue_y = (queue_Lon * (np.pi / 180) * (R * np.cos(Lat))) - (queue_Lon[0] * np.pi / 180 * (R * np.cos(Lat)))
plt.figure(9)
plt.plot(t,queue_x)
plt.title("N")
plt.xlabel('Cекунды')
plt.ylabel('N,км')
plt.grid(True)
plt.savefig('img/9.png')
plt.show()

plt.figure(10)
plt.plot(t,queue_y)
plt.title("E")
plt.xlabel('Cекунды')
plt.ylabel('E,км')
plt.grid(True)
plt.savefig('img/10.png')
plt.show()