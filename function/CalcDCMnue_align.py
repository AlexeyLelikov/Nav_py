import numpy as np

def CalcDCMnue_align(Pos : list,fb : np.array,om_ib : np.array) -> np.array:
    #nue->ned
    fb =np.array([fb[0],fb[2],-fb[1]])
    om_ib = np.array([om_ib[0],om_ib[2],-om_ib[1]])

    om_ie = 7.292115e-5

    cosB = np.cos(Pos[0])
    sinB = np.sin(Pos[0])
    cosL = np.cos(Pos[1])
    sinL = np.sin(Pos[1])

    Lat = Pos[0]

    (X,Y,Z) = LatLonAlt2XYZ(Pos[0],Pos[1],Pos[2])
    gravity = CalcGravity(X,Y,Z)

    gt = -cosL*sinB*gravity[0] - sinL*sinB*gravity[1] + cosB*gravity[2]
    gn = -(cosL*cosB*gravity[0] + sinL*cosB*gravity[1] + sinB*gravity[2])

    b1 = -fb
    b2 = om_ib
    b3 = cross(b1, b2)

    Mb = np.array([b1,b2,b3])

    Mn = np.array([[np.sin(Lat) / (gn * np.cos(Lat) + gt * np.sin(Lat)),gn / (gn * om_ie * np.cos(Lat) + gt * om_ie * np.sin(Lat)),0 ],
                  [0,0, 1 / (gn * om_ie * np.cos(Lat) + gt * om_ie * np.sin(Lat))],
                  [np.cos(Lat) / (gn * np.cos(Lat)+ gt * np.sin(Lat)),  -gt / (gn * om_ie * np.cos(Lat) + gt * om_ie * np.sin(Lat)), 0 ]] )

    R_nb = Mn @ Mb

    DCM = R_nb.T

    return DCM