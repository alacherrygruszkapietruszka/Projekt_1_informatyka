# Instrukcja obsługi programu skrypt.py
Program skrypt.py pozwala na przekształcanie współrzędnych geodezyjnych między różnymi układami odniesienia oraz formatami, wykorzystując dostępne elipsoidy geodezyjne. 

Wymagania:
Program dedykowany dla systemu operacyjnego windows. By program mógł działać konieczym jest posiadanie programu python 3.9 wraz z następującymi bibliotekami:
numpy (zainstalowane jako np): używane do operacji na macierzach i wektorach.
math: zawiera podstawowe funkcje matematyczne.
sys: zawiera metody i zmienne służące do modyfikowania wielu elementów środowiska wykonawczego języka Python

Elipsoidy geodezyjne:
Program obsługuje następujące elipsoidy geodezyjne:
WGS84: Elipsoida zdefiniowana w ramach systemu WGS84.
GRS80: Elipsoida zdefiniowana w ramach systemu GRS80.
mars: Elipsoida geodezyjna Marsa.

Dostępne transformacje:
XYZ_BLH: Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH.
BLH_XYZ: Przekształcenie współrzędnych BLH na współrzędne XYZ.
XYZ_NEU: Przekształcenie współrzędnych XYZ na współrzędne NEU.
BL_PL2000: Przekształcenie współrzędnych BLH na współrzędne PL-2000.
BL_PL1992: Przekształcenie współrzędnych BLH na współrzędne PL-1992.

Sposób użycia:
Uruchomienie programu wymaga podania parametrów przez wiersz poleceń.

Przykłady użycia:

Przekształcenie współrzędnych XYZ na współrzędne geodezyjne BLH:
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_BLH

Przekształcenie współrzędnych BLH na współrzędne XYZ:
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BLH_XYZ

Przekształcenie współrzędnych XYZ na współrzędne NEU:
python skrypt.py -plik dane.txt -elip WGS84 -funkcja XYZ_NEU

Przekształcenie współrzędnych BLH na współrzędne PL-2000:
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL2000

Przekształcenie współrzędnych BLH na współrzędne PL-1992:
python skrypt.py -plik dane.txt -elip WGS84 -funkcja BL_PL1992

Parametry:
-plik: nazwa pliku z danymi wejściowymi (należy podać pełną nazwę z rozszerzeniem, np. "dane.txt").
-elip: wybór elipsoidy (dozwolone wartości: "WGS84", "GRS80", "mars").
-funkcja: wybór rodzaju transformacji (dozwolone wartości: "XYZ_BLH", "BLH_XYZ", "XYZ_NEU", "BL_PL2000", "BL_PL1992").

Format danych wejściowych:
Dane wejściowe powinny być przechowywane w pliku tekstowym.
Współrzędne XYZ lub BLH powinny być oddzielone przecinkami.
Każda linia pliku powinna zawierać współrzędne dla jednego punktu w odpowiednim formacie.

Format danych wyjściowych:
Wyniki przekształcenia zostaną zapisane w pliku tekstowym o nazwie "WYNIK_NAZWA_FUNKCJI.txt", gdzie "NAZWA_FUNKCJI" to nazwa funkcji transformacji (np. "WYNIK_XYZ_BLH.txt").
Współrzędne będą oddzielone spacjami, a każda linia będzie zawierać współrzędne dla jednego punktu.
