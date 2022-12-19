import numpy as np
import pandas as pd
from scipy.io import loadmat
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

mat = loadmat('nav.binB_03.mat')
data_df = pd.DataFrame(mat['UAV'])
data = data_df.to_numpy().T
data = np.delete(data, 0, 1)
data = data[405500:,:]
Sensors = np.zeros((6,2),dtype = np.float64)
n = data.shape[0]
acc = data[:,9:12]
gyro = data[:,12:15]

GPSQuality = data[:,23]
GPS = data[:,17:23]
print(GPS.shape)
GPS[:,0:2] = GPS[:,0:2] * np.pi / 180

StatusMS = data[:, 33]
MagnSens = data[:, 30:33] * 10
MagnSens[:,1:3] = -1 * MagnSens[:,1:3]

dt = 1 / 100

print(GPS[1,0:2] * 180 / np.pi)



