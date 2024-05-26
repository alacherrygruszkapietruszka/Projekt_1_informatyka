# Instrukcja obsługi programu skrypt.py
Program skrypt.py pozwala na przekształcanie współrzędnych geodezyjnych między różnymi układami odniesienia oraz formatami, wykorzystując dostępne elipsoidy geodezyjne. 

# Wymagania:
Program dedykowany dla systemu operacyjnego windows. By program mógł działać konieczym jest posiadanie programu python 3.9 wraz z następującymi bibliotekami:\
numpy (zainstalowane jako np): używane do operacji na macierzach i wektorach.\
math: zawiera podstawowe funkcje matematyczne.\
sys: zawiera metody i zmienne służące do modyfikowania wielu elementów środowiska wykonawczego języka Python

# Elipsoidy geodezyjne:
Program obsługuje następujące elipsoidy geodezyjne:\
WGS84: Elipsoida zdefiniowana w ramach systemu WGS84.\
GRS80: Elipsoida zdefiniowana w ramach systemu GRS80.\
mars: Elipsoida geodezyjna Marsa.

# Dostępne transformacje:
XYZ_BLH: Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH.\
BLH_XYZ: Przekształcenie współrzędnych BLH na współrzędne XYZ.\
XYZ_NEU: Przekształcenie współrzędnych XYZ na współrzędne NEU.\
BL_PL2000: Przekształcenie współrzędnych BLH na współrzędne PL-2000.\
BL_PL1992: Przekształcenie współrzędnych BLH na współrzędne PL-1992.

# Sposób użycia:
Uruchomienie programu wymaga podania parametrów przez wiersz poleceń.

# Przykłady użycia:

Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH:\
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_BLH

Przekształcenie współrzędnych BLH na współrzędne XYZ:\
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BLH_XYZ

Przekształcenie współrzędnych XYZ na współrzędne NEU:\
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_NEU

Przekształcenie współrzędnych BLH na współrzędne PL-2000:\
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL2000

Przekształcenie współrzędnych BLH na współrzędne PL-1992:\
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL1992

# Parametry:
-plik: nazwa pliku z danymi wejściowymi (należy podać pełną nazwę z rozszerzeniem, np. "dane.txt").\
-elip: wybór elipsoidy (dozwolone wartości: "WGS84", "GRS80", "mars").\
-funkcja: wybór rodzaju transformacji (dozwolone wartości: "XYZ_BLH", "BLH_XYZ", "XYZ_NEU", "BL_PL2000", "BL_PL1992").

# Format danych wejściowych:
Dane wejściowe powinny być przechowywane w pliku tekstowym zapisanym w tym samym folderze w którym przechowywany jest program skrypt.py.\
Współrzędne XYZ lub BLH powinny być oddzielone przecinkami.\
Jednostka danych wyjścioweych jest zależna od wybranej funkcji:\
- XYZ_BLH - X[m], Y[m], Z[m]\
- BLH_XYZ - phi[stopnie dziesiętne], lambda[stopnie dziesiętne], H[m]\
- XYZ_NEU - Xa[m], Ya[m], Za[m], Xb[m], Yb[m], Zb[m]\
- BL_2000 - phi[stopnie dziesiętne], lambda[stopnie dziesiętne], m(stała - skala odwzorowania na południkach osiowych dla układu PL-2000)\
- BL_1992 - phi[stopnie dziesiętne], lambda[stopnie dziesiętne], m(stała - skala odwzorowania na południkach osiowych dla układu PL-1992)\
UWAGA! Korzystając z transformacji XYZ_NEU współrzędne środka oznaczone ...\
UWAGA! Program rozpoczyna wczytywanie danych  po pierwszych czterech linijkach nagłówka.\
Każda linia pliku (oprócz pierwszych czterech) powinna zawierać współrzędne dla jednego punktu w odpowiednim formacie.


Przykłady plkiów wejściowych:

Współrzedne geocentryczny ECEF stacji pemanentnej GNSS\
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu\
  X[m]         Y[m]        Z[m]\
#----------------------------------------------\
3664940.500,1409153.590,5009571.170\
3664940.510,1409153.580,5009571.167\
3664940.520,1409153.570,5009571.167\
3664940.530,1409153.560,5009571.168\
3664940.520,1409153.590,5009571.170\
3664940.514,1409153.584,5009571.166\
3664940.525,1409153.575,5009571.166\
3664940.533,1409153.564,5009571.169\
3664940.515,1409153.590,5009571.170\
3664940.514,1409153.584,5009571.169\
3664940.515,1409153.595,5009571.169\
3664940.513,1409153.584,5009571.171\

Współrzedne geodezyjne pomnika geodezji europejskiej\
...\
  phi[m]         lambda[m]        H[m]\
#----------------------------------------------\
52.243973,21.009177,109.91

Współrzedne geodezyjne ECEF stacji pemanentnej GNSS\
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu\
  phi[m]         lambda[m]         m\
#----------------------------------------------\

Współrzedne geodezyjne ECEF stacji pemanentnej GNSS\
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu\
  phi[m]         lambda[m]         m\
#----------------------------------------------\

Współrzedne geocentryczne ECEF stacji pemanentnej GNSS\
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu\
  Xa[m]         Ya[m]         Za[m]        Xb[m]         Yb[m]         Zb[m]\
#----------------------------------------------\



# Format danych wyjściowych:
Wyniki przekształcenia zostaną zapisane w pliku tekstowym o nazwie "WYNIK_NAZWA_FUNKCJI.txt", gdzie "NAZWA_FUNKCJI" to nazwa funkcji transformacji (np. "WYNIK_XYZ_BLH.txt").\
Współrzędne będą oddzielone spacjami, a każda linia będzie zawierać współrzędne dla jednego punktu.
Jednostka danych wyjściowych jest zależna od wybranej funkcji:\
- XYZ_BLH - phi[stopnie dziesiętne], lambda[stopnie dziesiętne], H[m]\
- BLH_XYZ - X[m], Y[m], Z[m]\
- XYZ_NEU - n[m], e[m], u[m]\
- BL_2000 - X2000[m], Y2000[m]\
- BL_1992 - X1992[m], Y1992[m]\


# Błędy:
- Program zwraca błąd przy podaniu niewłaściwej liczby argumentów i prosi o podanie wszystkich,
- Program zwraca błąd przy podaniu parametrów bez wartości,
- Program zwraca błąd w przypadku podania nieprawidłowego modelu elipsoidy,
- Program zwraca błąd w przypadku podania nieobsługiwanej transformacji i prosi o podanie jednej z możliwych opcji,
- Program zwraca błąd kiedy plik wejściowy nie może być znaleziony i prosi o podanie innego pliku, wprowadzenie jego lokalizacji lub sprawdzenie nazwy pliku,
- Program rozpoczyna wczytywanie danych  po pierwszych czterech linijkach nagłówka.
- Program zwraca błąd w przypadku podania pliku o niewłaściwym formacie (musi byc txt)


