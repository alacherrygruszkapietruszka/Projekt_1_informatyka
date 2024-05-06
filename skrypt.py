from math import sin, cos, sqrt, atan, atan2, degrees, radians
import numpy as np
o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        def xyz2flh(self, X, Y, Z):
            flh = []
            for X,Y,Z in zip(X,Y,Z):
                p = np.sqrt(X**2 + Y**2)
                fi = np.arctan(Z / (p * (1 - self.e2)))
                while True:
                    N = self.Npu(fi)
                    h = p / np.cos(fi) - N
                    fip = fi     #fip - fi poprzednie, fi - fi nowe
                    fi = np.arctan(Z / (p * (1 - N * self.e2 / (N + h))))
                    if abs(fip - fi) < (0.000001/206265):
                        break
            
                lam = np.arctan2(Y, X)
                flh.extend([np.rad2deg(fi), np.rad2deg(lam), h])
            return(flh)