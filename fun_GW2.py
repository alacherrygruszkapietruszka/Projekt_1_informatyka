# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 17:51:09 2023

@author: weron
"""

from math import *
import numpy as np

# PARAMETRY ELIPSOIDY GRS80
a_grs= 6378137     
e2_grs = 0.00669438002290 
#%%  MOJE DANE
Xa = 3667475.000 - (9*10)
Ya = 1100280.000 + (9*20)
Za = 5084025.000 + (9*30)
s_AB = 30000.000 + (9*1000.000)
A_AB = radians(140)
alfa_AB = radians(95.0000 + 9*5) 
z_AB = radians(90)    

#%%  HIRVONEN ----- DEF 
# promien krzywizny w 1 werykale
def Np(f,a,e2):
    N = a/np.sqrt(1-e2*np.sin(f)**2)
    return(N)

def hirvonen(X,Y,Z,a,e2):
    '''
    Parameters: X, Y, Z, a, e2 (pkt A)
    ----------
    Returns: f, l, h (pkt A)
    -------
    '''
    #promień równoleżnika
    p = np.sqrt(X**2 + Y**2)
    f = np.arctan(Z / (p*(1-e2)))
    while True: 
        N = Np(f,a,e2)
        h = (p/np.cos(f)) - N
        fs = f
        f = np.arctan(Z/(p * (1 - (e2 * (N / (N + h))))))
        if np.abs(fs-f) < (0.000001/206265):
            break
    l = np.arctan2(Y,X) 
    return(f,l,h)

def dms(x):
    '''
    Przeliczanie wartosci z radianow do miary czasowej
    ---------
    Parameters: x_rad
    ---------
    Returns: x_stopnie_minuty_sekundy
    '''
    sig = ' '
    if x<0:
        sig = '-'
        x = abs(x)
    x = x * 180 / pi
    d = int(x)
    m = int(60 * (x - d))
    s = (x - d - m/60) * 3600
    print(f'{sig}{d:3d}{chr(176)}{abs(m):2d}\'{abs(s):7.5f}\"')
    
# promien przekroju normalnego podluznego
def Mp(a,e2,f):
    M = (a * (1 - e2))/(np.sqrt(1-e2*np.sin(f)**2)**3)
    return(M)

# sredni promien krzywizny
def Rp(M,N):
    R = np.sqrt(M*N)
    return(R)

#kontrola
def flh2XYZ(h,f,l,a,e2):
    '''
    Przeliczanie wspolrzednych prostokatnych na wspolrzedne geodezyjne
    ----------
    Parameters: h, f, l, a, e2 
    ----------
    Returns: X, Y, Z

    '''
    N = a/np.sqrt(1-e2*np.sin(f)**2)
    X = (N+h)*np.cos(f)*np.cos(l)
    Y = (N+h)*np.cos(f)*np.sin(l)
    Z = (N*(1-e2)+h)*np.sin(f)
    return(X, Y, Z)

#%%  KIVIOJ ----- DEF 
def kivioj(f,l,A,s,a,e2):
    '''
    Parameters: f, l, Az, s, a, e2 (pkt A)
    ----------
    Returns: f, l, Az  (pkt B)
    -------
    '''
    n=int(s/1000)                         #n - liczba kawałków, s/1000 - cała długosc przez dokładnosc
    ds= s/n                               #delta s
    for i in range(n):  
        M = Mp(a,e2,f)                    #główne promienie krzyeizny w pkt wyjsciowym
        N = Np(f,a,e2)                    # to samo
        df = ds*np.cos(A)/M               #przyblizenie przyrostow szerokosci
        dA = (ds*np.sin(A)*np.tan(f))/N   #-:- azymutu
        fm = f+df/2                       #srednia wartosc szerokosci i azymutu w pkt srodkowym m
        Am=A+dA/2
        Mm= Mp(a,e2,fm)
        Nm = Np(fm,a,e2)
        df = ds*np.cos(Am)/Mm             #ostateczne przyrosty 
        dl = (np.sin(Am)*ds)/(Nm*np.cos(fm))
        dA =(np.sin(Am)*np.tan(fm)*ds)/Nm 
        f = f+df                          #wspolrzedne i azymuut pkt koncowego
        l = l+dl
        A = A+dA 
    A = A + pi 
    if A > 2* pi:          
        A = A - 2*pi      
    return(f,l,A)

#%%  ZAMIANA MIARY CZASOWEJ NA RADIANY
def z_cz_na_r(stp,mi,sek):
    '''
    Parameters: wartosc_stp_mi_sek
    ----------
    Returns: wartosc_rad
    ----------
    '''
    w_s = stp + mi/60 + sek/3600
    w_r = radians(w_s)
    return(w_r)

#wynik_w_r = z_cz_na_r(stp,mi,sek)
#podstawic wartosci za stopnie min, sekundy
#%%  VINCENTY ----- DEF 
def vincenty(fa,la,fb,lb,a,e2):
    '''
    Parameters: fa, la, fb, lb, a, e2
    ----------
    Returns s_AB, A_AB, A_BA
    -------
    '''
    #krótsza półos elipsoidy
    b = a*np.sqrt(1-e2)
    #splaszczenie elipsoidy
    fl = 1 - b/a    
    #roznica  dlugosci geodezyjnej
    dL = lb - la
    #szerokosc zredukowana
    Ua = np.arctan((1-fl)*np.tan(fa))
    Ub = np.arctan((1-fl)*np.tan(fb))
    L = dL                          
    while True:
        sin_sigma = np.sqrt((np.cos(Ub)*np.sin(L))**2 + (np.cos(Ua)*np.sin(Ub)-np.sin(Ua)*np.cos(Ub)*np.cos(L))**2)
        cos_sigma = np.sin(Ua)*np.sin(Ub) + np.cos(Ua)*np.cos(Ub)*np.cos(L)
        #odleglosc katowa miedzy punktami na sferze
        sigma = np.arctan(sin_sigma/cos_sigma)
        # alfa - azymut linii geodezyjnej na rowniku
        sin_alfa = (np.cos(Ua)*np.cos(Ub)*np.sin(L))/sin_sigma  
        cos2_alfa = 1 - sin_alfa**2
        #sigma m - odleglosc katowa na sferze od rownika do pkt srodkowego linii geodezyjnej
        cos_2_sigma_m = cos_sigma - ((2*np.sin(Ua)*np.sin(Ub))/cos2_alfa)
        C = fl/16*cos2_alfa*(4+fl*(4-3*cos2_alfa))
        Ls = L
        #roznica dlugosci na sferze pomocnoczej
        L = dL + (1-C)*fl*sin_alfa*(sigma+C*sin_sigma*(cos_2_sigma_m + C*cos_sigma*(-1+2*((cos_2_sigma_m)**2))))
        if abs(L-Ls)<(0.000001/206265):
            break
    u2 = ((a**2 - b**2)/b**2)*cos2_alfa
    A = 1 + (u2/16384)*(4096+u2*(-768+u2*(320-175*u2)))
    B = (u2/1024)*(256+u2*(-128+u2*(74-47*u2)))
    dsigma = B*sin_sigma*(cos_2_sigma_m + (1/4)*B*(cos_sigma*(-1+2*((cos_2_sigma_m)**2))-(1/6)*B*cos_2_sigma_m*(-3+4*((sin_sigma)**2))*(-3+4*((cos_2_sigma_m)**2))))
    #dlugosc linii geodezyjnej
    s = b*A*(sigma-dsigma)
    AAB = np.arctan2(np.cos(Ub)*np.sin(L) , np.cos(Ua)*np.sin(Ub)-np.sin(Ua)*np.cos(Ub)*np.cos(L))   
    ABA = np.arctan2(np.cos(Ua)*np.sin(L) , (-np.sin(Ua)*np.cos(Ub)+np.cos(Ua)*np.sin(Ub)*np.cos(L)) ) + pi
    return(s,AAB, ABA) 

#%%  NEU ----- DEF
def saz2neu(s,alfa,z):   #mam s, alfa, z: chce x,y,z
    dx = np.array([s * np.sin(z) * np.cos(alfa),    #dn
                   s * np.sin(z) * np.sin(alfa),    #de
                   s * np.cos(z) ])                 #du
    return(dx)

def R_w_neu(f,l):
    R = np.array([[-np.sin(f) * np.cos(l), -np.sin(l), np.cos(f) * np.cos(l)],
                  [-np.sin(f) * np.sin(l), np.cos(l), np.cos(f) * np.sin(l)],
                  [np.cos(f), 0, np.sin(f)]])
    return(R)

def neu2XYZ(dx,f,l):
    R = R_w_neu(f,l)
    dX = R @ dx
    return(dX)

def XYZ2neu(dx, f, l):
    R = R_w_neu(f,l)
    dX = R.T @ dx
    return(dX)

def neu2saz(dX):
    dn, de, du = dX[0], dX[1], dX[2]
    s = np.sqrt(dn**2 + de**2 + du**2)
    alfa = np.arctan(de/dn)
    z = np.arccos(du/s) 
    return(s,alfa,z)

#%%  GAUSS KRUGER
# dlugosc luku poludnika
def sigma(f, a, e2):
    A_0 = 1 - (e2/4) - ((3*(e2**2))/64) - ((5*(e2**3))/256)
    A_2 = (3/8) * (e2 + (e2**2)/4 + (15*(e2**3))/128)
    A_4 = (15/256) * (e2**2 + (3*(e2**3))/4)
    A_6 = (35*(e2**3))/3072
    sigma = a * (A_0 * f - A_2 * np.sin(2*f) + A_4 * np.sin(4*f) - A_6 * np.sin(6*f))
    return (sigma)

def gauss_kruger(a, e2, f, l, l0):
    '''
    Parameters: a, e2, f, l, l0
    ----------
    Returns: xgk, ygk
    -------
    '''
    b2 = a**2 * (1-e2)
    e_prim2 = (a**2 - b2)/b2
    delta_l = l - l0
    t = np.tan(f)
    eta2 = e_prim2 * ((np.cos(f))**2) 
    N = Np(f,a,e2)
    si = sigma(f, a, e2)
    xgk = si + ((delta_l)**2)/2 * N * np.sin(f) * np.cos(f) * (1 + ((delta_l)**2)/12 * (np.cos(f))**2 * (5 - t**2 + 9 * eta2 + 4 * (eta2)**2) + ((delta_l)**4)/360 * (np.cos(f))**4 * (61 - 58*t**2 + t**4 + 270*eta2 - 330*eta2*t**2))
    ygk = delta_l * N * np.cos(f) * (1 + ((delta_l)**2)/6 * (np.cos(f))**2 * (1 - t**2 + eta2) + ((delta_l)**4)/120 * (np.cos(f))**4 * (5 - 18 * t**2 + t**4 + 14 * eta2 - 58 * eta2 * t**2))
    return(xgk, ygk)

def gk2pl92(xgk, ygk, m92 = 0.9993):
    x92 = xgk * m92 - 5300000
    y92 = ygk * m92 + 500000
    return(x92, y92)

def gk2pl2000(xgk, ygk, nr_strefy, m2000 = 0.999923):
    x2000 = xgk * m2000 
    y2000 = ygk * m2000 + 500000 + nr_strefy * 1000000
    return(x2000, y2000)

def pl20002gk(x_2000, y_2000, nr_strefy, m2000 = 0.999923):
     x_gk = x_2000/m2000
     y_gk = (y_2000 - 500000 - nr_strefy * 1000000)/m2000
     return (x_gk, y_gk)
 
def pl922gk(x_92, y_92, m92 = 0.9993):
    xgk = (x_92 + 5300000) / m92
    ygk = (y_92 - 500000) / m92
    return (xgk, ygk)

def f1(xgk, a, e2):
    A_0 = 1 - (e2/4) - ((3*(e2**2))/64) - ((5*(e2**3))/256)
    f1 = xgk/(a*A_0)
    while True:
        fs = f1
        si = sigma(f1,a,e2)
        f1 = f1 + (xgk - si)/(a * A_0)
        if np.abs(f1 - fs) < (0.000001 / 206265):
            break
    return (f1)

def gauss_kruger_odw(xgk,ygk,a,e2, l0):
    '''
    Parameters: xgk, ygk, a, e2, l0
    ----------
    Returns: f, l
    -------
    '''
    f_1 = f1(xgk, a, e2)
    N_1 = Np(f_1,a,e2)
    t_1 = np.tan(f_1)
    M_1 = Mp(a,e2,f_1)
    b2 = a**2 * (1-e2)
    e_prim2 = (a**2 - b2)/b2
    eta2_1 =  e_prim2 * (np.cos(f_1))**2
    fi = f_1 - (((ygk)**2) * t_1)/(2*M_1*N_1) * (1 - ((ygk)**2)/(12 * N_1**2)*(5 + 3 * (t_1)**2 + eta2_1 - 9 * eta2_1 * (t_1**2) - 4 * (eta2_1)**2) + ((ygk)**4)/(360 * (N_1)**4) * (61 + 90 * (t_1)**2 + 45 * (t_1)**4))
    l = l0 + ygk/(N_1 * np.cos(f_1)) * (1 - ygk**2/(6*(N_1)**2) * (1 + 2 * t_1**2 + eta2_1) + ygk**4 /(120 * N_1**4) * (5 + 28 * t_1**2 + 24 * t_1**4 + 6 * eta2_1 + 8 * eta2_1 * (t_1**2)))
    return(fi, l)

#%%  REDUKCJE ----- DEF

# skala odzworowania
def mgk(xgk, ygk, a, e2): 
    f = f1(xgk, a, e2)
    R = np.sqrt(Np(f,a,e2) * Mp(a,e2,f))
    m = 1 + (ygk**2)/(2 * R**2) + ygk**4/(24 * R**4)
    return (m)

# skala układu
def m_uklad(mgk, m0_uklad):
    m_uklad = mgk * m0_uklad
    return(m_uklad)

# redukcja dlugosci skosnej na elipsoide
def r_ds_na_elipsoide(fa, fb, a, e2, H_A, H_B, s_pom):
    fm = (fa + fb)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    delta_H_AB = H_B - H_A
    s_0 = np.sqrt(((s_pom**2) - (delta_H_AB**2)) / ((1 + H_A/Rm) * (1 + H_B/Rm))) 
    s_elip = 2 * Rm * np.arcsin(s_0/(2*Rm))
    return(s_elip)

# redukcja odleglosci
def red_d(fa, fb, a, e2, s_elip, ygk_A, ygk_B):
    fm = (fa + fb)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    r = s_elip * ((ygk_A**2 + ygk_A * ygk_B + ygk_B**2)/(6*(Rm**2)))
    sgk = s_elip + r
    return(sgk)

# wsp srednie
def wsp_srednie(xa,xb,ya,yb):
    x_sr = (xa + xb) / 2
    y_sr = (ya + yb) / 2
    return(x_sr, y_sr)
#%% REDUKCJE 2 ----- DEF

# zbieznosc poludnikow z fi i lambda
def zbiez_pol_fl(l, l0, f, a, e2):
    b2 = a**2 * (1-e2)
    e_prim2 = (a**2 - b2)/b2
    delta_l = l - l0
    eta2 = e_prim2 * ((np.cos(f))**2)
    t = np.tan(f)
    gamma = delta_l * np.sin(f) + delta_l**3/3 * np.sin(f) * (np.cos(f))**2 * (1 + 3*eta2 + 2*eta2**2) + delta_l**5 / 15 * np.sin(f) * np.cos(f)**4 * (2 - t*2)
    return(gamma)

# zbieznosc poludnikow z x i y
def zbiez_pol_xy(xgk, ygk, a, e2):
    f_1 = f1(xgk, a, e2)
    N_1 = Np(f_1,a,e2)
    t_1 = np.tan(f_1) 
    b2 = a**2 * (1-e2)
    e_prim2 = (a**2 - b2)/b2
    eta2_1 =  e_prim2 * (np.cos(f_1))**2
    gamma = ygk/N_1 * t_1 * (1 - ygk**2/(3*N_1**2) * (1 + t_1**2 - eta2_1 - 2*eta2_1**2) + ygk**4/(15*N_1**4) * (2 + 5*t_1**2 + 3*t_1**4))
    return(gamma)

# redukcja kierunku              
def red_kier_AB(xagk, xbgk, yagk, ybgk, a, e2, l0):
    fi_a, l_a = gauss_kruger_odw(xagk, yagk, a, e2, l0)
    fi_b, l_b = gauss_kruger_odw(xbgk, ybgk, a, e2, l0)
    fm = (fi_a + fi_b)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    delta_AB = ((xbgk - xagk)*(2*yagk + ybgk))/(6*Rm**2)
    return(delta_AB)
    
def red_kier_BA(xagk, xbgk, yagk, ybgk, a, e2, l0):
    fi_a, l_a = gauss_kruger_odw(xagk, yagk, a, e2, l0)
    fi_b, l_b = gauss_kruger_odw(xbgk, ybgk, a, e2, l0)
    fm = (fi_a + fi_b)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    delta_BA = ((xagk - xbgk)*(2*ybgk + yagk))/(6*Rm**2)
    return(delta_BA)

def red_kier(xagk, xbgk, yagk, ybgk, a, e2, l0):
    fi_a, l_a = gauss_kruger_odw(xagk, yagk, a, e2, l0)
    fi_b, l_b = gauss_kruger_odw(xbgk, ybgk, a, e2, l0)
    fm = (fi_a + fi_b)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    delta_AB = ((xbgk - xagk)*(2*yagk + ybgk))/(6*Rm**2)
    delta_BA = ((xagk - xbgk)*(2*ybgk + yagk))/(6*Rm**2)
    return(delta_AB, delta_BA)

def alfa1(xa,ya,xb,yb):
    dx = xb - xa
    dy = yb - ya
    alfa_ab = atan2(dy,dx)
    alfa_ba = alfa_ab +  pi
    return(alfa_ab,alfa_ba)
# redukcja azymutu
def red_az(xagk, xbgk, yagk, ybgk, a, e2, l0):
    alfa_AB = atan2((ybgk - yagk),(xbgk - xagk))
    alfa_BA = alfa_AB + pi
    gamma_A = zbiez_pol_xy(xagk, yagk, a, e2)
    gamma_B = zbiez_pol_xy(xbgk, ybgk, a, e2)
    delta_AB = red_kier_AB(xagk, xbgk, yagk, ybgk, a, e2, l0)
    delta_BA = red_kier_BA(xagk, xbgk, yagk, ybgk, a, e2, l0)
    A_AB = alfa_AB + gamma_A + delta_AB
    A_BA = alfa_BA + gamma_B + delta_BA
    return(A_AB, A_BA)

# redukcja długosci AB
def red_d_xy(xgk_A, xgk_B, ygk_A, ygk_B, l0, xau, xbu, yau, ybu, m0, a, e2):
    fi_a, l_a = gauss_kruger_odw(xgk_A, ygk_A, a, e2, l0)
    fi_b, l_b = gauss_kruger_odw(xgk_B, ygk_B, a, e2, l0)
    fm = (fi_a + fi_b)/2
    Mm = Mp(a,e2,fm)
    Nm = Np(fm, a, e2)
    Rm = Rp(Mm,Nm)
    s = np.sqrt((xbu - xau)**2 + (ybu - yau)**2)
    sgk = s/m0
    r = sgk * ((ygk_A**2 + ygk_A * ygk_B + ygk_B**2)/(6*(Rm**2)))
    sAB = sgk - r
    return(r, sAB)

#%% TRANSFORMACJA ----- DEF
def wek_A(r_prim):
    '''
    Parameters: wspolrzedne przestrzenne punktow w ukladzie pierwotnym
    ----------
    Funkcja uzywa wspolrzednych w ukladzie pierwotnym r',
    nastepnie zwraca jedna macierz A, ktora ma ilosc elementow rowna ilosci
    punktow *3. Kolejne 3 wiersze odpowiadaja jednemu punktowi

    Returns: macierz A dla wszystkich punktow
    -------
    '''
    x_prim = r_prim[:, 0]
    y_prim = r_prim[:, 1]
    z_prim = r_prim[:, 2]
    A_ = []
    for i in range(0, r_prim.shape[0]):
        A = np.array([[x_prim[i], 0, -z_prim[i], y_prim[i], 1, 0, 0],
                      [y_prim[i], z_prim[i], 0, -x_prim[i], 0, 1, 0],
                      [z_prim[i], -y_prim[i], x_prim[i], 0, 0, 0, 1]])
        A_.extend(A)
    return (np.array(A_))


def wek_L(r_prim, r_bis):
    '''
    Parameters: wspolrzedne w ukladzie pierwotnym r' i wtornym r''
    ----------
    Funkcja tworzy macierz L dla r' i r''. Ilosc wierszy macierzy L
    zalezy od ilosci punktow w ukladzie r' i r'' (ilosc wierszy = ilosc punktow*3)
    
    Returns: macierz L dla wszystkich punktow
    -------
    '''
    x_prim = r_prim[:, 0]
    y_prim = r_prim[:, 1]
    z_prim = r_prim[:, 2] 
    x_bis = r_bis[:, 0]
    y_bis = r_bis[:, 1]
    z_bis = r_bis[:, 2]
    L_ = []
    for i in range(0, r_prim.shape[0]):
        L = np.array([[x_bis[i] - x_prim[i]],
                      [y_bis[i] - y_prim[i]],
                      [z_bis[i] - z_prim[i]]])
        L_.extend(L)
    return(np.array(L_))

def trans(alfa, beta, gamma, kappa, xyzgrs80, r0):
    '''
    Obliczenie wspolrzednych punktu polozonego na elipsoidzie Krasowskiego.
    
    Parameters:
    alfa : parametr transformacji (kat)
    beta : parametr transformacji (kat)
    gamma : parametr transformacji (kat)
    kappa : [niemianowane]
    xyzgrs80 : wspolrzedne XYZ punktu na elipsoidzie GRS'80
    r0 : wektor translacji r0 = [X0, Y0, Z0]

    Returns:
    XYZ_pkt : wspolrzedne prostokatne szukanego punktu 

    '''
    d_B = np.array([[0, gamma, -beta],
                  [-gamma, 0, alfa],
                  [beta, -alfa, 0]])
    d_M = np.diag([kappa, kappa, kappa])
    E = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]])
    XYZ_pkt = (E + d_M + d_B) @ xyzgrs80 + r0
    return(XYZ_pkt)

def transf(x, y, z, kx, ky, kz, alfa, beta, gamma, x0, y0, z0):
    r1 = np.array([x, y, z])
    M = np.array([[kx, gamma, -beta],
                  [-gamma, ky, alfa],
                  [beta, -alfa, kz]])
    r0 = np.array([x0, y0, z0])
    r2 = r1 + M @ r1 + r0
    return(r2)