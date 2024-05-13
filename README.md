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
UWAGA! Program rozpoczyna wczytywanie danych  po pierwszych czterech linijkach nagłówka.\
Każda linia pliku (oprócz pierwszych czterech) powinna zawierać współrzędne dla jednego punktu w odpowiednim formacie.\

- Przykład plkiu wejściowego:

 Współrzedne geocentryczny ECEF stacji pemanentnej GNSS
Obserwatorium Astronomiczno-Geodezyjne w Józefosławiu
  X[m]         Y[m]        Z[m]
# ----------------------------------------------------
3664940.500,1409153.590,5009571.170
3664940.510,1409153.580,5009571.167
3664940.520,1409153.570,5009571.167
3664940.530,1409153.560,5009571.168
3664940.520,1409153.590,5009571.170
3664940.514,1409153.584,5009571.166
3664940.525,1409153.575,5009571.166
3664940.533,1409153.564,5009571.169
3664940.515,1409153.590,5009571.170
3664940.514,1409153.584,5009571.169
3664940.515,1409153.595,5009571.169
3664940.513,1409153.584,5009571.171


# Format danych wyjściowych:
Wyniki przekształcenia zostaną zapisane w pliku tekstowym o nazwie "WYNIK_NAZWA_FUNKCJI.txt", gdzie "NAZWA_FUNKCJI" to nazwa funkcji transformacji (np. "WYNIK_XYZ_BLH.txt").\
Współrzędne będą oddzielone spacjami, a każda linia będzie zawierać współrzędne dla jednego punktu.

# Błędy:
- Program zwraca błąd przy podaniu niewłaściwej liczby argumentów i prosi o podanie wszystkich,\
- Program zwraca błąd przy podaniu parametrów bez wartości,\
- Program zwraca błąd w przypadku podania nieprawidłowego modelu elipsoidy,\
- Program zwraca błąd w przypadku podania nieobsługiwanej transformacji i prosi o podanie jednej z możliwych opcji,\
- Program zwraca błąd kiedy plik wejściowy nie może być znaleziony i prosi o podanie innego pliku, wprowadzenie jego lokalizacji lub sprawdzenie nazwy pliku,\
- Program rozpoczyna wczytywanie danych  po pierwszych czterech linijkach nagłówka.
- Program zwraca błąd w przypadku podania pliku o niewłaściwym formacie (musi byc txt)
gi

