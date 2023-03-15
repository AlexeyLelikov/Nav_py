import numpy as np
def LatLonAlt2XYZ(Lat : float,Lon : float,Alt : float) -> tuple:

    a_e = 6378137
    ecc = 0.081819190842622
    ksi = 1.0 / np.sqrt(1 - ecc * ecc * np.sin(Lat) * np.sin(Lat))
    X = (a_e * ksi + Alt) * np.cos(Lat) * np.cos(Lon)
    Y = (a_e * ksi + Alt) * np.cos(Lat) * np.sin(Lon)
    Z = (a_e * ksi * (1 - ecc * ecc) + Alt) * np.sin(Lat)
    return X, Y, Z