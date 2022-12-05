import numpy as np
from SkewSymmMatr import SkewSymmMatr
def cross(a : np.array ,b : np.array) -> np.array:
    axb = SkewSymmMatr(a) @ b