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

        def Npu(self, fi):     #promien krzywizny w I wertykale
            N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
            return(N)

        def sigma(self, fi):
            A0 = 1 - (self.e2/4) - (3*(self.e2)**2)/64 -  (5*(self.e2)**3)/256
            A2 = 3/8 * (self.e2 + (self.e2)**2/4 + 15*(self.e2)**3/128)
            A4 = 15/256 * ( (self.e2)**2 + (3*((self.e2)**3))/4 )
            A6 = 35 * (self.e2)**3 / 3072
            sigma = self.a * ( A0 * fi - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi) )
            return(sigma)        

        def xyz2flh(self, X, Y, Z, output = 'dec_degree'):
            """
            Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
            na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
            W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
            Parameters
            ----------
            X, Y, Z : FLOAT
                 współrzędne w układzie orto-kartezjańskim, 

            Returns
            -------
            lat
                [stopnie dziesiętne] - szerokość geodezyjna
            lon
                [stopnie dziesiętne] - długośc geodezyjna.
            h : TYPE
                [metry] - wysokość elipsoidalna
            output [STR] - optional, defoulf 
                dec_degree - decimal degree
                dms - degree, minutes, sec
            """
            r   = sqrt(X**2 + Y**2)           # promień
            lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
            lat = 0
            while abs(lat_prev - lat) > 0.000001/206265:    
                lat_prev = lat
                N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
                h = r / cos(lat_prev) - N
                lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
            lon = atan(Y/X)
            N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
            h = r / cos(lat) - N       
            if output == "dec_degree":
                return degrees(lat), degrees(lon), h 
            elif output == "dms":
                lat = self.deg2dms(degrees(lat))
                lon = self.deg2dms(degrees(lon))
                return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
            else:
                raise NotImplementedError(f"{output} - output format not defined")
                

            def plh2xyz(self, phi, lam, h):
                """
                Algorytm zamiany współrzędnych geodezyjnych: długość, szerokość i wysokośc elipsoidalna (phi, lam, h) na 
                współrzędne ortokartezjańskie (x, y, z).
                Parameters
                ----------
                phi
                    [stopnie dziesiętne] - szerokość geodezyjna
                lam
                    [stopnie dziesiętne] - długośc geodezyjna.
                h : TYPE
                    [metry] - wysokość elipsoidalna

                Returns
                -------
                X, Y, Z : FLOAT
                     [metry] - współrzędne w układzie orto-kartezjańskim, 
                
                """
                phi = radians(phi)
                lam = radians(lam)
                Rn = self.a/sqrt(1 - self.ecc2 * sin(phi)**2)
                q = Rn * self.ecc2 * sin(phi)
                x = (Rn + h) * cos(phi) * cos(lam)
                y = (Rn + h) * cos(phi) * sin(lam)
                z = (Rn + h) * sin(phi) - q 
                return x, y, z
            
            
            def xyz2neu(self, x, y, z, x_0, y_0, z_0):
                """
                Algorytm przeniesiena współrzędnych ortokaretzjanskich (x, y, z) na współrzędne w układzie
                horyzontalnum neu.

                Parameters
                ----------
                x, y, z : FLOAT
                    [metry] - współrzędne w układzie orto-kartezjańskim
                x_0, y_0, z_0 : FLOAT
                    [metry] - współrzędne w układzie orto-kartezjańskim
                    
                Returns
                -------
                N, E, U : FLOAT
                    [metry] - współrzędne w układzie horyzontalnym

                """
                phi, lam, _ = [radians(coord) for coord in self.xyz2plh(x_0, y_0, z_0)]
                
                R = np.array([[-sin(lam), -sin(phi)*cos(lam), cos(phi)*cos(lam)],
                              [ cos(lam), -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                              [        0,           cos(phi),          sin(phi)]])
                
                xyz_t = np.array([[x - x_0],
                                  [y - y_0],
                                  [z - z_0]])
                [[E], [N], [U]] = R.T @ xyz_t
                
                return N, E, U

           
            
            def PL2000(self, f, l, m=0.999923):
                """
                Algorytm zamiany współrzędnych geodezyjnych na wspórzędne w układzie współrzędnych
                PL-2000.

                Parameters
                ----------
                f : FLOAT
                    [stopnie dziesiętne] - szerokosć geodezyjna
                l : FLOAT
                    [stopnie dziesiętne] - długosć geodezyjna
                m: FLOAT
                    stała - skala odwzorowania na południkach osiowych dla układu PL-2000
                Returns
                -------
                result : LIST
                    lista zawierająca kolejne współrzędne x, y w układzie PL-2000

                """
                result = []
                for f, l in zip(f,l):
                   l0 = 0 
                   strefa = 0
                   if l >=np.deg2rad(13.5) and l <= np.deg2rad(16.5):
                       strefa = 5
                       la0 = np.deg2rad(15)
                   elif l >np.deg2rad(16.5) and l <= np.deg2rad(19.5):
                       strefa = 6
                       l0 = np.deg2rad(18)
                   elif l >np.deg2rad(19.5) and l <= np.deg2rad(22.5):
                       strefa =7
                       l0 = np.deg2rad(21)
                   elif l >np.deg2rad(22.5) and l <= np.deg2rad(25.5):
                       strefa = 8
                       l0 = np.deg2rad(24)
                   b2 = (self.a**2) * (1-self.ep2)   #krotsza polos
                   e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
                   dl = l - l0
                   t = np.tan(f)
                   ni = np.sqrt(e2p * (np.cos(f))**2)
                   N = self.Np(f)
                   sigma = self.sigma(f)
                   XGK20 = sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1 + ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
                   YGK20 = (dl*N* np.cos(f)) * (1+(((dl)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                   X2000 = XGK20 * m 
                   Y2000 = YGK20 * m + strefa*1000000 + 500000
                   result.append([X2000, Y2000])
               
                   return result

            def PL1992(self, f, l, m = 0.9993):
                """
                Algorytm zamiany współrzędnych geodezyjnych na wspórzędne w układzie współrzędnych
                PL-1992.

                Parameters
                ----------
                f : FLOAT
                    [stopnie dziesiętne] - szerokosć geodezyjna
                l : FLOAT
                    [stopnie dziesiętne] - długosć geodezyjna
                m: FLOAT
                    stała - skala odwzorowania na południkach osiowych dla układu PL-1992
                Returns
                -------
                result : LIST
                    lista zawierająca kolejne współrzędne x, y w układzie PL-1992

                """
                result = []
                lam0 = (19*np.pi)/180
                for f, l in zip(f,l):
                    b2 = (self.a**2) * (1-self.ep2)   #krotsza polos
                    e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
                    dlam = l - lam0
                    t = np.tan(f)
                    ni = np.sqrt(e2p * (np.cos(f))**2)
                    N = self.Np(f)
            
                    sigma = self.sigma(f)
            
                    xgk = sigma + ((dlam**2)/2)*N*np.sin(f)*np.cos(f) * ( 1+ ((dlam**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
                    ygk = (dlam*N* np.cos(f)) * (1+(((dlam)**2/6)*(np.cos(f))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                    
                    x92 = xgk*m - 5300000
                    y92 = ygk*m + 500000
                       
                    result.append([x92, y92])
                   
                    return result 
        
        if __name__ == "__main__":
            # utworzenie obiektu
            geo = Transformacje(model = "wgs84")
            # dane XYZ geocentryczne
            X = 3664940.500; Y = 1409153.590; Z = 5009571.170
            phi, lam, h = geo.xyz2plh(X, Y, Z)
            print(phi, lam, h)
            # phi, lam, h = geo.xyz2plh2(X, Y, Z)
            # print(phi, lam, h)
                    