import numpy as np
from Savage import Savage
import matplotlib.pyplot as plt
from LatLonAlt2XYZ import LatLonAlt2XYZ
from CalcGravity import CalcGravity
from SkewSymmMatr import SkewSymmMatr
from cross import cross
from Calc_F_C_N import Calc_F_C_N
from Calc_C_N_E import Calc_C_N_E
from twoPiBound import twoPiBound
from LatLonAlt2XYZ import LatLonAlt2XYZ
from DCM_bn import DCM_bn
from CalcDCMnue_align import CalcDCMnue_align
from EulerFromDCMnue import EulerFromDCMnue


# Чтение файла
def find_znach(buf_str, flag):
    if buf_str != '':
        # разбиваем строку
        split_str = buf_str.split('  ')
        # print(split_str)

        time = split_str[0]

        # три акселерометра
        accel_1 = split_str[7]
        accel_2 = split_str[8]
        accel_3 = split_str[9]
        # гироскопы
        sav_1 = split_str[10]
        sav_2 = split_str[11]
        sav_3 = split_str[12]

        GPS_lat = split_str[18]
        GPS_lon = split_str[19]
        p_fromFile = split_str[1]  # Тангаж
        r_fromFile = split_str[2]  # Крен
        h_fromFile = split_str[3]  # Курс
        Ve = split_str[13]
        Vn = split_str[14]
        # print(accel_1,accel_2,accel_3, sav_1, sav_2,sav_3)

        flag = True
        return [time, accel_1, accel_2, accel_3, sav_1, sav_2, sav_3, GPS_lat, GPS_lon, p_fromFile, r_fromFile,
                h_fromFile, Ve, Vn, flag]

    else:
        flag = False
        return ['' for i in range(9)] + [flag]


f = open('CN5_20130907_Mi8_Vertical_T_Kazan-Tumen.txt', 'r')
# f=open('/content/drive/MyDrive/Alg/05-ConvArRT-BINS-dt-1_ar05c6u1_1.txt', 'r')
# f=open('/content/drive/MyDrive/Alg/20210428200455_binout1_nku4.txt', 'r')

# полезные измерения, начинаются с 3669-й строки  (до этого по GPS координатам и измерениям нули)

deg2rad = np.pi / 180
Lat = []
Lon = []

Lat_fromFile = []
Lon_fromFile = []

roll = []
pitch = []
heading = []

a_x = [] # E
a_y = [] # N
a_z = [] # u

w_x = [] # E
w_y = [] # N
w_z = [] # U

time_ar = []
p_fromFile = []
r_fromFile = []
h_fromFile = []

Ve_fromFile = [] # E
Vn_fromFile = [] # N
Lat_fromFile.append(55.6)
Lon_fromFile.append(49.2)
ft = 0.01
flag_EOF = True  # устанавливаем флаг в истину, т.к. мы только что открыли файл и находимся в его начале
x = 0
it = 0
tbuf = 0
cl = 0
j = 0
time_fromFile = []

# flag_EOF = False
for i in range(3670):  # 36 сек
    buf_str = f.readline()
    t, accel_1, accel_2, accel_3, sav_1, sav_2, sav_3, GPS_lat, GPS_lon, p_fromFile_str, r_fromFile_str, h_fromFile_str, Ve, Vn, flag_EOF = find_znach(
        buf_str, flag_EOF)
    time_fromFile.append(ft * j / 60)
    Ve_fromFile.append(float(Ve))
    Vn_fromFile.append(float(Vn))

    if float(GPS_lat) != 0:
        Lat_fromFile.append(float(GPS_lat))
    else:
        Lat_fromFile.append(Lat_fromFile[i - 1])
    if float(GPS_lon) != 0:
        Lon_fromFile.append(float(GPS_lon))
    else:
        Lon_fromFile.append(Lon_fromFile[i - 1])
    p_fromFile.append(float(p_fromFile_str))  # Тангаж
    r_fromFile.append(float(r_fromFile_str))  # Крен
    h_fromFile.append(float(h_fromFile_str))  # Курс
    j = j + 1

while (flag_EOF) and (cl != 60 * 60 / 0.01):
    buf_str = f.readline()  # читаем строку из файла, чтобы потом определить значения градуса, минут и секунд
    # print(buf_str)
    t, accel_1, accel_2, accel_3, sav_1, sav_2, sav_3, GPS_lat, GPS_lon, p_fromFile_str, r_fromFile_str, h_fromFile_str, Ve, Vn, flag_EOF = find_znach(
        buf_str, flag_EOF)

    if flag_EOF:

        time_ar.append(float(t))
        if float(GPS_lat) != 0:
            Lat_fromFile.append(float(GPS_lat))
        else:
            Lat_fromFile.append(Lat_fromFile[i - 1])
        if float(GPS_lon) != 0:
            Lon_fromFile.append(float(GPS_lon))
        else:
            Lon_fromFile.append(Lon_fromFile[i - 1])
        a_x.append(float(accel_1))
        a_y.append(float(accel_2))
        a_z.append(float(accel_3))
        w_x.append(float(sav_1) * deg2rad)
        w_y.append(float(sav_2) * deg2rad)
        w_z.append(float(sav_3) * deg2rad)
        p_fromFile.append(float(p_fromFile_str))  # Тангаж
        r_fromFile.append(float(r_fromFile_str))  # Крен
        h_fromFile.append(float(h_fromFile_str))  # Курс
        time_fromFile.append(ft * j / 60)
        Ve_fromFile.append(float(Ve))
        Vn_fromFile.append(float(Vn))
        j = j + 1
        cl = cl + 1
f.close()

Lat_fromFile = Lat_fromFile[1:]
Lon_fromFile = Lon_fromFile[1:]
N = 363670
T = N / 100 # c
t = np.linspace(0,T,N)


# ================================================== Данные =========================================================
dt = 1 / 100
(a_x , a_y, a_z) = (a_y , a_z, a_x) # ENU -> NUE
(w_x , w_y, w_z) = (w_y , w_z, w_x) # ENU -> NUE
acc = np.hstack((np.array(a_x).reshape(-1,1) , np.array(a_y).reshape(-1,1), np.array(a_z).reshape(-1,1)))  * dt # NUE
gyro = np.hstack((np.array(w_x).reshape(-1,1) , np.array(w_y).reshape(-1,1), np.array(w_z).reshape(-1,1))) * dt # NUE
print(acc)
n = acc.shape[0]
print(n)
s = 10000
n = n - s
print(n)

# ============================================= Выставка =========================================================
Lat = Lat_fromFile[s] * np.pi / 180
Lon = Lon_fromFile[s] * np.pi / 180
Alt = 0
Lat_old = Lat
Alt_old = Alt
W_NUE = np.array([Ve_fromFile[s],Vn_fromFile[s],0]) # ENU
W_NUE_old = W_NUE
Roll = r_fromFile[s] * np.pi / 180
Pitch = p_fromFile[s] * np.pi / 180
Heading = h_fromFile[s] * np.pi / 180
C_B_N = DCM_bn(Heading, Pitch, Roll)

NavState = np.array([Lat,Lon,Alt,          # Широта, долгота, геодезическая высота
                   W_NUE[0],W_NUE[1],W_NUE[2],               # Скорость на север, вверх, на восток
                   Roll,Pitch,Heading,   # Крен, тангаж, курс
                   0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,           # Ошибки масштабных коэффициентов акселерометров, смещения нуля акселерометров, ошибки масштабных коэффициентов гироскопов
                   Lat_old,Alt_old,              # Предыдущие широта и высота
                   W_NUE_old[0],W_NUE_old[1],W_NUE_old[2]])

Sensors = np.zeros((6,2),dtype = np.float64)
queue_Lat = np.zeros(n)
queue_Lon = np.zeros(n)
queue_pitch = np.zeros(n)
queue_roll = np.zeros(n)
queue_heading = np.zeros(n)
queue_Alt = np.zeros(n)
queue_W_N = np.zeros(n)
queue_W_U = np.zeros(n)
queue_W_E = np.zeros(n)

## ======================= Работа ===========================================================================
for i in range(s,n):
    Sensors[0:3,1] = acc[i,:] # выходной сигнал акселерометра по осям чувствительности (x,y,z) на текущем такте расчетов
    Sensors[0:3,0] = acc[i-1,:]
    Sensors[3:6,0] = gyro[i-1,:] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на предыдущем такте расчетов
    Sensors[3:6,1] = gyro[i,:] # выходной сигнал гироскопов (скорость) по оси (x,y,z) на текущем такте расчетов
    NavState[2] = 0
    NavState[4] = 0
    queue_Lat[i-s] = NavState[0]
    queue_Lon[i-s] = NavState[1]
    queue_roll[i-s] = NavState[6]
    queue_pitch[i-s] = NavState[7]
    queue_heading[i - s] = NavState[8]
    queue_W_N[i-s] = NavState[3]
    queue_W_E[i-s] = NavState[5]
    queue_Alt[i-s] = NavState[2]
    (NavState, delta_v, delta_alpha, gt, C_B_N) = Savage(NavState, Sensors, dt, C_B_N)

queue_Lat = queue_Lat * (180 / np.pi)
queue_Lon = queue_Lon * (180 / np.pi)
queue_pitch = queue_pitch * (180 / np.pi)
queue_roll = queue_roll * (180 / np.pi)
queue_heading = queue_heading * (180 / np.pi)

n = queue_Lat.shape[0]
N = n
T = N / 100 # c
t = np.linspace(0,T,N)


plt.figure(1)
plt.plot(t,queue_Lat)
plt.plot(t,Lat_fromFile[s:-3670],'r')
plt.title("Широта")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(2)
plt.plot(t,queue_Lon)
plt.plot(t,Lon_fromFile[s:-3670],'r')
plt.title("Долгота")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(3)
plt.plot(t,queue_pitch)
plt.plot(t,p_fromFile[s:-3670],'r')
plt.title("Тангаж")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(4)
plt.plot(t,queue_roll)
plt.plot(t,r_fromFile[s:-3670],'r')
plt.title("Крен")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(5)
plt.plot(t,queue_heading)
plt.plot(t,h_fromFile[s:-3670],'r')
plt.title("Курс")
plt.xlabel('Cекунды')
plt.ylabel('град')
plt.grid(True)
plt.show()

plt.figure(6)
plt.plot(t,queue_W_E)
plt.title("Скорость E")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()

plt.figure(7)
plt.plot(t,queue_W_N)
plt.title("Скорость N")
plt.xlabel('Cекунды')
plt.ylabel('м/c')
plt.grid(True)
plt.show()