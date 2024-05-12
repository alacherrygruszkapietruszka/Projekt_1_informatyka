from math import sin, cos, sqrt, atan, atan2, degrees, radians
import numpy as np
import argparse as args, argparse
from argparse import ArgumentParser
import os

o = object()

class Transformacje:
    def __init__(self, elipsoida):
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
        self.a = elipsoida[0]
        self.e2 = elipsoida[1]
        """
        if model == "WGS84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "GRS80":
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
        """
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

    def xyz2flh(self, X, Y, Z, output='dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : ARRAY
             współrzędne w układzie orto-kartezjańskim, 
    
        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : ARRAY
            [metry] - wysokość elipsoidalna
        output [STR] - optional, default 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        blh = []
        r = np.sqrt(X**2 + Y**2)  # promień
        lat_prev = np.arctan(Z / (r * (1 - self.e2)))  # pierwsze przybliilizenie
        lat = np.zeros_like(lat_prev)
        while np.any(np.abs(lat_prev - lat) > 0.000001 / 206265):
            lat_prev = lat
            N = self.a / np.sqrt(1 - self.e2 * np.sin(lat_prev)**2)
            h = r / np.cos(lat_prev) - N
            lat = np.arctan((Z / r) * (((1 - self.e2 * N / (N + h))**(-1))))
        lon = np.arctan2(Y, X)
        N = self.a / np.sqrt(1 - self.e2 * (np.sin(lat))**2)
        h = r / np.cos(lat) - N
        if output == "dec_degree":
            return np.degrees(lat), np.degrees(lon), h
        elif output == "dms":
            lat = self.deg2dms(np.degrees(lat))
            lon = self.deg2dms(np.degrees(lon))
            blh = blh.extend([lat, lon, h])
            return blh
        else:
            raise NotImplementedError(f"{output} - output format not defined")

                

    def flh2xyz(self, phi, lam, h):
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
        Rn = self.a/sqrt(1 - self.e2 * sin(phi)**2)
        q = Rn * self.e2 * sin(phi)
        x = (Rn + h) * cos(phi) * cos(lam)
        y = (Rn + h) * cos(phi) * sin(lam)
        z = (Rn + h) * sin(phi) - q 
        return x, y, z
            
    def get_dXYZ(self, xa, ya, za, xb, yb, zb):
        '''
        funkcja liczy macierz różnicy współrzednych punktów A i B, która jest potrzebna do obliczenia macierzy neu

        Parametry
        ----------
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        dXYZ : ARRAY
            macierz różnicy współrzędnych

        '''
        dXYZ = np.array([xb-xa, yb-ya, zb-za])
        return(dXYZ)
    
    
    def rneu(self, f, l):
        '''
        Funkcja tworzy macierz obrotu R, która jest potrzebna do obliczenia macierzy neu

        Parametry
        ----------
        f : FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna..
        l : FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.

        Returns
        -------
        R ARRAY
            macierz obrotu R
             
        '''
        f=np.radians(f)
        l=np.radians(l)
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                      [-np.sin(f)*np.sin(l),  np.cos(l), np.cos(f)*np.sin(l)],
                      [np.cos(f),             0,         np.sin(f)          ]])
        return(R)        
    def xyz2neu(self, f, l, xa, ya, za, xb, yb, zb):
        '''
        Układ współrzędnych horyzontalnych – układ współrzędnych astronomicznych, w którym oś główną stanowi 
        lokalny kierunek pionu, a płaszczyzną podstawową jest płaszczyzna horyzontu astronomicznego. 
        Biegunami układu są zenit i nadir. Ich położenie na sferze niebieskiej zależy od współrzędnych geograficznych 
        obserwatora oraz momentu obserwacji, tak więc współrzędne horyzontalne opisują jedynie chwilowe położenie ciała niebieskiego.

        Parametry
        ----------
        f : FLOAT
            [stopnie dziesiętne] - szerokość geodezyjna..
        l : FLOAT
            [stopnie dziesiętne] - długośc geodezyjna.
        XA, YA, ZA, XB, YB, ZB: FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        n , e, u : STR
            współrzędne horyzontalne
            

        '''
        neu = []
        dX = Transformacje.get_dXYZ(self, xa, ya, za, xb, yb, zb)
        R = Transformacje.rneu(self, f,l)
        neu = R.T @ dX
        n = neu[0];   e = neu[1];   u = neu[2]
        n = "%.16f"%n; e = "%.16f"%e; u="%.16f"%u
        dlugosc = []
        xx = len(n); dlugosc.append(xx)
        yy = len(e); dlugosc.append(yy)
        zz = len(u); dlugosc.append(zz)
        P = 50
        
        while xx < P:
            n = str(" ") + n
            xx += 1
        
        while yy < P:
            e = str(" ") + e
            yy += 1
            
        while zz < P:
            u = str(" ") + u
            zz +=1
        neu.extend([n, e, u])
        return(neu)
    
           
            
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
           b2 = (self.a**2) * (1-self.e2)   #krotsza polos
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
            b2 = (self.a**2) * (1-self.e2)   #krotsza polos
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
     
    def open_and_save(self, plik, funkcja):

        with open(plik, 'r') as file:
            lines = file.readlines()
    
        # Pominięcie nagłówka składającego się zarówno z liter, jak i cyfr
        data_start_index = 0
        for i, line in enumerate(lines):
            if not any(c.isalpha() for c in line.strip()) or all(c.isdigit() or c in '.,' for c in line.strip()):
                data_start_index = i
                break
    
        data = np.genfromtxt(lines[data_start_index:], delimiter=",")
    
        np.set_printoptions(precision=5, suppress=True)
    
        if funkcja == "XYZ_BLH":
            X = data[:,0]
            Y = data[:,1]
            Z = data[:,2]
            blh = self.xyz2flh(X, Y, Z)
            header = "Transformacja XYZ -> BLH"
            # Zapisujemy dane po trzy w wierszu
            blh_str = '\n'.join(','.join(f'{x:.5f}' for x in row) for row in blh)
            with open(os.path.join(os.getcwd(), f"WYNIK_{funkcja.upper()}.txt"), 'w') as f:
                f.write(header + '\n')
                f.write(blh_str)
        elif funkcja == "BLH_XYZ":
            phi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            h = data[:,2]
            XYZ = self.flh2xyz(phi, lam, h)
            header = "Transformacja BLH -> XYZ"
            np.savetxt(os.path.join(os.getcwd(), f"WYNIK_{funkcja.upper()}.txt"), XYZ, delimiter=",", header=header, fmt='%.5f')
        elif funkcja == "BL_PL1992":
            f1 = np.deg2rad(data[:,0])
            l1 = np.deg2rad(data[:,1])
            result92 = self.PL1992(f1, l1)
            header = "Transformacja BLH -> PL1992"
            np.savetxt(os.path.join(os.getcwd(), f"WYNIK_{funkcja.upper()}.txt"), result92, delimiter=",", header=header, fmt='%.5f')
        elif funkcja == "BL_PL2000":
            f = np.deg2rad(data[:,0])
            l = np.deg2rad(data[:,1])
            result00 = self.PL2000(f, l)
            header = "Transformacja BLH -> PL2000"
            np.savetxt(os.path.join(os.getcwd(), f"WYNIK_{funkcja.upper()}.txt"), result00, delimiter=",", header=header, fmt='%.5f')
        elif funkcja == "XYZ_NEU":
            f = data[0,0]
            l = data[0,1]
            X0 = data[0,2]
            Y0 = data[0,3]
            Z0 = data[0,5]
            X = data[0,6]
            Y = data[0,7]
            Z = data[0,8]
            neu = self.xyz2neu(f, l, X, Y, Z, X0, Y0, Z0)
            header = "Transformacja XYZ -> NEU"
            np.savetxt(os.path.join(os.getcwd(), f"WYNIK_{funkcja.upper()}.txt"), neu, delimiter=",", header=header, fmt='%.5f')

 
           
if __name__ == "__main__":
    
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-plik", type = str, help = "Podaj nazwę pliku, w którym znajdują się jako dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument('-elip', type=str, help="Podaj jedną z wskazanych elipsoid: GRS80, WGS84, mars")
        parser.add_argument("-funkcja", type = str, help = "Wybierz transformację jaką chcesz obliczyć: 'XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU' ")
        args = parser.parse_args()
    except SyntaxError:
        print(f"Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")
           
    elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'mars':[3396900.000, 0.008725863396028901]}
    funkcja = {'XYZ_BLH':'xyz2flh', 'BLH_XYZ': 'flh2xyz', 'XYZ_NEU': 'xyz2neu', 'BL_PL2000': 'PL2000', 'BL_PL1992': 'PL1992'}
    
    try:
        geo = Transformacje(elip[args.elip.upper()])
        finito = geo.open_and_save(args.plik, args.funkcja.upper())
        print("Zapisano")
    except KeyError:
        print(f"Podana funkcja/elipsoida nie istnieją, proszę upewnij się, że korzystasz z istniejących elipsoid")
    except AttributeError:
        print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartosci.")
    except FileNotFoundError:
        print("Nie znaleziono takiego pliku. Proszę spróbować wprowadzić inny plik.")
    except IndexError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
    except ValueError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
                 
