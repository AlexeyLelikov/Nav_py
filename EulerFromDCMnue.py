import numpy as np

def threeaxisrot(r11, r12, r21, r31, r32):
    r1 = np.arctan2(r11,r12)
    r2 = np.arcsin(r21)
    r3 = np.arctan2(r31,r32)
    return r1,r2,r3

def EulerFromDCMnue(dcm):
    (r1, r2, r3) = threeaxisrot(dcm[0, 1],dcm[0,0], -dcm[0, 2],dcm[1, 2], dcm[2, 2])
    return (r1, r2, r3)
