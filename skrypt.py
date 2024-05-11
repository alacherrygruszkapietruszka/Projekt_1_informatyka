from math import sin, cos, sqrt, atan, atan2, degrees, radians
import numpy as np
import argparse as args, argparse
from argparse import ArgumentParser


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
            
        return(n, e, u)
    
           
            
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
     
    def file_open(self, name):
        """
        Wczytanie pliku .txt i wyodrębnienie podanych w nim współrzędnych
        za pomocą pętli for. Odczytane dane dodajemy do list.
        
        Parameters:
        -----------
        name : str
            Nazwa pliku do odczytania.
    
        Returns:
        --------
        X : list[float]
            Lista zawierająca współrzędne X.
        Y : list[float]
            Lista zawierająca współrzędne Y.
        Z : list[float]
            Lista zawierająca współrzędne Z.
        """
        X = []
        Y = []
        Z = []
        with open(name, 'r') as file:
            tab = np.genfromtxt(file, delimiter = "," , dtype ='<U20', skip_header = 4)
            for line in tab:
                x = line[0]
                X.append(float(x[0]))
                Y.append(float(x[1]))
                Z.append(float(x[2]))
        return X, Y, Z

    def file_save92(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> PL-1992 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        X92 = []
        Y92 = []
        for a, b, c in zip(X, Y, Z):
            x92, y92 = Transformacje.PL1992(a, b, c)           
            X92.append(x92)
            Y92.append(y92)
            
        plik=open(file_out,"w")
        plik.write(f'# PL-1992---------------------------------------------- \n')
        plik.write(f'  X[m]         Y[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b in zip(X92,Y92):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            plik.write(f'{a},   {b} \n')
        plik.close()
    
    def file_save00(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> PL-2000 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        X00 = []
        Y00 = []
        for a, b, c in zip(X, Y, Z):
            x00, y00 = Transformacje.PL2000(a, b, c)
            X00.append(x00)
            Y00.append(y00)
            
        plik=open(file_out,"w")
        plik.write(f'# PL-2000---------------------------------------------- \n')
        plik.write(f'  X[m]         Y[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b in zip(X00,Y00):
            a = f'{a:7.3f}'
            b = f'{b:7.3f}'
            plik.write(f'{a},   {b} \n')
        plik.close()

    
    def file_saveFLH(self, X, Y, Z, file_out):
        """
        Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji XYZ -> BLH 
        do pliku wyjsciowego file_out w formacie .txt
        """
        
        F = []
        L = []
        H = []
        for a, b, c in zip(X, Y, Z):
            f, l, h = Transformacje.xyz2flh(a, b, c)
            F.append(degrees(f))
            L.append(degrees(l))
            H.append(h)
            
        plik=open(file_out,"w")
        plik.write(f'  B[d]         L[d]         H[m] \n')
        plik.write(f'# ----------------------------------------------------- \n')
        for a,b,c in zip(F,L,H):
            a = f'{a:7.4f}'
            b = f'{b:7.4f}'
            c = f'{c:7.3f}'
            plik.write(f'{a},      {b},      {c} \n')
        plik.close()
            
    def file_saveNEU(self, X, Y, Z, file_out):
        """
        Funkcja pozwala na zapisanie wyników transformacji XYZ -> NEU do pliku.
    
        Parameters:
        -----------
        X : float
            Współrzędna X w układzie ortokartezjańskim.
        Y : float
            Współrzędna Y w układzie ortokartezjańskim.
        Z : float
            Współrzędna Z w układzie ortokartezjańskim.
        file_out : str
            Nazwa pliku do zapisania.
    
        Returns:
        --------
        None
        """
        dX = self.get_dXYZ(X[0], Y[0], Z[0], X[1], Y[1], Z[1])
        R = self.rneu(self.f, self.l)
        neu = R.T @ dX
        n = "%.16f" % neu[0]
        e = "%.16f" % neu[1]
        u = "%.16f" % neu[2]
        dlugosc = []
        xx = len(n)
        dlugosc.append(xx)
        yy = len(e)
        dlugosc.append(yy)
        zz = len(u)
        dlugosc.append(zz)
        P = 50
    
        while xx < P:
            n = str(" ") + n
            xx += 1
    
        while yy < P:
            e = str(" ") + e
            yy += 1
    
        while zz < P:
            u = str(" ") + u
            zz += 1
    
        with open(file_out, "w") as plik:
            plik.write(f'  N[m]         E[m]         U[m] \n')
            plik.write(f'# ----------------------------------------------------- \n')
            plik.write(f'{n},   {e},      {u} \n')
            
        def file_saveXYZ(self, B, L, H, file_out):
            """
            Funkcja ta pozwala zapisać współrzędne otrzymane po transformacji BLH -> XYZ 
            do pliku wyjsciowego file_out w formacie .txt
            """
            
            X = []
            Y = []
            Z = []
            for a, b, c in zip(B, L, H):
                x, y, z = Transformacje.plh2xyz(a, b, c)
                X.append(x)
                Y.append(y)
                Z.append(z)
                
            plik=open(file_out,"w")
            plik.write(f'  X[m]         Y[m]         Z[m] \n')
            plik.write(f'# ----------------------------------------------------- \n')
            for a,b,c in zip(X,Y,Z):
                a = f'{a:7.3f}'
                b = f'{b:7.3f}'
                c = f'{c:7.3f}'
                plik.write(f'{a},      {b},      {c} \n')
            plik.close()
 
           
if __name__ == "__main__":
    geo = Transformacje("grs80")
    
    
    parser = ArgumentParser()
    parser = argparse.ArgumentParser(description="Podaj plik")
    parser.add_argument("-plik", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
    parser.add_argument('-elip', '--elip', type=str, help="Podaj jedną z wskazanych elipsoid: GRS80, WGS84, mars")
    parser.add_argument('-neu', '--neu', type=str, help="Podaj nazwe pliku wynikiowego dla neu z rozszerzeniem txt")
    parser.add_argument('-xa', '--xa', type=float)
    parser.add_argument('-ya', '--ya', type=float)
    parser.add_argument('-za', '--za', type=float)
    parser.add_argument('-xb', '--xb', type=float)
    parser.add_argument('-yb', '--yb', type=float)
    parser.add_argument('-zb', '--zb', type=float)
    args = parser.parse_args()
 
    elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'mars':[3396900.000, 3376097.80585952]}
 
    geo = Transformacje(model = args.elip)
    f, l, h = geo.xyz2flh(args.xa, args.ya, args.za)
    n, e, u = geo.xyz2neu(f, l, args.xa, args.ya, args.za, args.xb, args.yb, args.zb)
     
     
    geo.file_saveNEU(args.neu, n, e, u)
     
    n = float(n)
    e = float(e)
    u = float(u)
     
     
    try:
        geo = Transformacje(elip[args.elip.upper()])
        finito = geo.plik(args.plik, args.funkcja.upper())
        print("Zapisano")
    except KeyError():
        print(f"Podana funkcja/elipsoida nie istnieją, proszę upewnij się, że korzystasz z istniejących elipsoid")
    except AttributeError:
        print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartosci.")
    except FileNotFoundError:
        print("Nie znaleziono takiego pliku. Proszę spróbować wprowadzić inny plik.")
    except IndexError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
    except ValueError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
                 
