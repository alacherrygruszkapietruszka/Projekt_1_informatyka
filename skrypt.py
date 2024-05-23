from math import sin, cos, sqrt, atan, atan2, degrees, radians
import numpy as np
import sys

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

    def Npu(self, fi):
        """
        Promień krzywizny w I wertykale
        Parameters
        ----------
        fi: FLOAT
            szerokosc geocentrycza [radiany]
        
        Returns
        -------
        N FLOAT
            promień krzywizny w I wertykale [metry]
        """
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
        return N

    def sigma(self, fi):
        A0 = 1 - (self.e2/4) - (3*(self.e2)**2)/64 -  (5*(self.e2)**3)/256
        A2 = 3/8 * (self.e2 + (self.e2)**2/4 + 15*(self.e2)**3/128)
        A4 = 15/256 * ( (self.e2)**2 + (3*((self.e2)**3))/4 )
        A6 = 35 * (self.e2)**3 / 3072
        sigma = self.a * ( A0 * fi - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi) )
        return sigma        

    def xyz2flh(self, X, Y, Z):
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
        lat STR
            [stopnie dziesiętne] - szerokość geodezyjna
        lon STR
            [stopnie dziesiętne] - długośc geodezyjna.
        h : STR
            [metry] - wysokość elipsoidalna
        """
        r = np.sqrt(X**2 + Y**2)
        lat_prev = np.arctan(Z / (r * (1 - self.e2)))
        lat = np.zeros_like(lat_prev)
        while np.any(np.abs(lat_prev - lat) > 0.000001 / 206265):
            lat_prev = lat
            N = self.a / np.sqrt(1 - self.e2 * np.sin(lat_prev)**2)
            h = r / np.cos(lat_prev) - N
            lat = np.arctan((Z / r) * (((1 - self.e2 * N / (N + h))**(-1))))
        lon = np.arctan2(Y, X)
        N = self.a / np.sqrt(1 - self.e2 * (np.sin(lat))**2)
        h = r / np.cos(lat) - N
        
        lat = np.rad2deg(lat)
        lon = np.rad2deg(lon)
        return f'{lat} {lon} {h} \n'
        
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
       X, Y, Z : STR
            [metry] - współrzędne w układzie orto-kartezjańskim
       """

        phi = np.radians(phi)
        lam = np.radians(lam)
        Rn = self.a / np.sqrt(1 - self.e2 * np.sin(phi) ** 2)
        q = Rn * self.e2 * np.sin(phi)
        x = (Rn + h) * np.cos(phi) * np.cos(lam)
        y = (Rn + h) * np.cos(phi) * np.sin(lam)
        z = (Rn * (1 - self.e2) + h) * np.sin(phi)
        return f'{x} {y} {z} \n'

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
        return dXYZ
    
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

        f = np.radians(f)
        l = np.radians(l)
        R = np.array([[-np.sin(f)*np.cos(l), -np.sin(l), np.cos(f)*np.cos(l)],
                      [-np.sin(f)*np.sin(l),  np.cos(l), np.cos(f)*np.sin(l)],
                      [np.cos(f),             0,         np.sin(f)]])
        return R

    def xyz2neu(self, xa, ya, za, xb, yb, zb):
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

        p = np.sqrt(X0**2+Y0**2)
        f = np.arctan(Z0/(p*(1-self.e2)))
        l = np.arctan(Y0/X0)
        dX = self.get_dXYZ(xa, ya, za, xb, yb, zb)
        R = self.rneu(f, l)
        neu = R.T @ dX
        n = neu[0];   e = neu[1];   u = neu[2]
        n = "%.16f"%n; e = "%.16f"%e; u="%.16f"%u
        return f'{str(neu)} \n'
    
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
        X2000, Y2000 STR
            f'string zawierający współrzędne x i y w układzie PL-2000

        """

        l0 = 0 
        strefa = 0
        if l >= radians(13.5) and l <= radians(16.5):
            strefa = 5
            l0 = radians(15)
        elif l > radians(16.5) and l <= radians(19.5):
            strefa = 6
            l0 = radians(18)
        elif l > radians(19.5) and l <= radians(22.5):
            strefa =7
            l0 = radians(21)
        elif l > radians(22.5) and l <= radians(25.5):
            strefa = 8
            l0 = radians(24)
        b2 = (self.a**2) * (1-self.e2)
        e2p = ( self.a**2 - b2 ) / b2
        l = radians(l)
        dl = l - l0
        t = np.tan(f)
        ni = e2p * (np.cos(f))**2
        N = self.Npu(f)
        sigma = self.sigma(f)
        XGK20 = sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1 + ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*ni + 4*(ni**2))  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*ni - 330*ni*(t**2)))
        YGK20 = (dl*N* np.cos(f)) * (1+((((dl)**2)/6)*(np.cos(f))**2) *(1-(t**2)+ni)+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*ni-58*ni*(t**2)) )
        X2000 = XGK20 * m 
        Y2000 = YGK20 * m + strefa*1000000 + 500000
        return f'{X2000} {Y2000} \n'

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
           stała - skala odwzorowania na południkach osiowych dla układu PL-2000
       Returns
       -------
       X92, Y92 STR
           f'string zawierający współrzędne x i y w układzie PL-1992

       """
        l = radians(l)
        lam0 = radians(19)
        b2 = (self.a**2) * (1-self.e2)
        e2p = ( self.a**2 - b2 ) / b2
        dl = l - lam0
        t = np.tan(f)
        ni = np.sqrt(e2p * (np.cos(f))**2)
        N = self.Npu(f)
        sigma = self.sigma(f)
        XGK92 = sigma + ((dl**2)/2)*N*np.sin(f)*np.cos(f) * ( 1 + ((dl**2)/12)*(np.cos(f))**2 * ( 5 - (t**2)+9*ni + 4*(ni**2))  + ((dl**4)/360)*(np.cos(f)**4) * (61-58*(t**2)+(t**4) + 270*ni - 330*ni*(t**2)))
        YGK92 = (dl*N* np.cos(f)) * (1+((((dl)**2)/6)*(np.cos(f))**2) *(1-(t**2)+ni)+((dl**4)/120)*(np.cos(f)**4)*(5-18*(t**2)+(t**4)+14*ni-58*ni*(t**2)) )
        X92 = XGK92*m - 5300000
        Y92 = YGK92*m + 500000
        return f'{X92} {Y92} \n'

if __name__ == "__main__":

    elipsoidy = {
        'WGS84':[6378137.000, 0.00669438002290],
        'GRS80':[6378137.000, 0.00669438002290],
        'mars':[3396900.000, 0.008725863396028901]}

    funkcja = {'XYZ_BLH':'xyz2flh',
               'BLH_XYZ': 'flh2xyz',
               'XYZ_NEU': 'xyz2neu',
               'BL_PL2000': 'PL2000',
               'BL_PL1992': 'PL1992'}

    argumenty = sys.argv[1:]
    
    if '-plik' not in argumenty or '-elip' not in argumenty or '-funkcja' not in argumenty:
        raise Exception('Nie podano wszystkich argumentów: -plik, -elip, -funkcja. Podaj wszytkie argumenty.')
    
    try:
        elip = argumenty[argumenty.index('-elip') + 1]
        trans_wsp = argumenty[argumenty.index('-funkcja') + 1]
        plik = argumenty[argumenty.index('-plik') + 1]
    except IndexError:
        raise Exception('Nie podano wartosci dla wszytkich wymaganych parametrów')
        
    try:
        elipsoida = elipsoidy[elip]
        geo = Transformacje(elipsoida)
    except KeyError:
        raise Exception('Podano niewłasciwy typ modelu elipsoidy. Podaj jeden z dostępnych: GRS80, WGS84, mars.')

    transformacje = ['XYZ_BLH', 'BLH_XYZ', 'BL_PL2000', 'BL_PL1992', 'XYZ_NEU']
    
    if trans_wsp not in funkcja:
        raise Exception('Skrypt nie obsluguje podanej transformacji. Podaj jedną z możliwych: XYZ_BLH:xyz2flh,BLH_XYZ:flh2xyz,XYZ_NEU:xyz2neu,BL_PL2000: PL2000,BL_PL1992: PL1992')
        
    try:
        with open(plik, 'r') as f, open(f"WYNIK_{trans_wsp.upper()}.txt", 'w') as wynik:
            linie = f.readlines()
            linie = linie[4:]
            for index, linia in enumerate(linie): 
                linia = linia.strip()
                if trans_wsp in ['XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU']:
                    x_str, y_str, z_str = linia.split(',')
                    x, y, z = float(x_str), float(y_str), float(z_str)
                    if trans_wsp == 'XYZ_BLH':
                        wynik.write(geo.xyz2flh(x, y, z))
                    elif trans_wsp == 'BLH_XYZ':
                        wynik.write(geo.flh2xyz(x, y, z))
                    elif trans_wsp == 'XYZ_NEU':
                        if index == 0:
                            X0, Y0, Z0 = x, y, z
                            continue
                        neu = geo.xyz2neu(X0, Y0, Z0, x, y, z)
                        wynik.write(' '.join(neu) + '\n')
                else:
                    fi_str, lam_str = linia.split(',')
                    fi, lam = float(fi_str), float(lam_str)
                    if trans_wsp == 'BL_PL2000':
                        wynik.write(geo.PL2000(fi, lam))
                    elif trans_wsp == 'BL_PL1992':
                        wynik.write(geo.PL1992(fi, lam))
    except FileNotFoundError:
        raise Exception('Podany plik nie istnieje. Podaj inny plik, sprawdz jego lokalizacje lub sprawdz nazwę podanego pliku.')
    except (KeyError, IndexError, ValueError):
        raise Exception('Nieodpowiedni format pliku.')
    except AttributeError:
        raise print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartosci.")
    
    print('Zapisano. Wyniki znajdują się w pliku WYNIK_<funkcja>.txt')
