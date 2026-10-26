# Lab 2 – ciało odbijające się od powierzchni w polu grawitacyjnym

Ten plik wyjaśnia, jak działa program [lab2.py](lab2.py) i skąd biorą się wyniki.

## 1. Na czym polega zadanie

Ciało (punkt materialny o masie 1 kg) zostaje wyrzucone z punktu startowego z prędkością początkową.
Spada w jednorodnym polu grawitacyjnym ($g = 10\ \text{m/s}^2$) i odbija się od powierzchni $z = f(x, y)$.
W przykładzie jest to paraboloida $z = x^2 + y^2$, czyli „miska”.

Trzeba:

1. zasymulować pierwsze **5 odbić** bez strat energii, podać punkty uderzenia i czasy między nimi
   oraz narysować powierzchnię z trajektoriami,
2. policzyć w punktach uderzenia energię kinetyczną, potencjalną i całkowitą i sprawdzić, czy energia jest stała,
3. powtórzyć to samo, gdy przy każdym odbiciu ciało traci część energii (współczynnik $k = 0.8$).

Dane z przykładu (na górze pliku `lab2.py`):

| wielkość | wartość |
|---|---|
| start $(x_0, y_0, z_0)$ | $(-5,\ 0,\ 200)$ |
| prędkość $(v_{0x}, v_{0y}, v_{0z})$ | $(2,\ 4,\ 0)$ |
| powierzchnia | $f(x, y) = x^2 + y^2$ |
| straty | $k = 0.8$ |

## 2. Fizyka w skrócie

### Lot między odbiciami

Jedyną siłą jest grawitacja, więc przyspieszenie wynosi $[0, 0, -g]$.
W poziomie ciało porusza się ze stałą prędkością, a w pionie jak w rzucie pionowym:

$$
\begin{aligned}
x(t) &= x_0 + v_{0x} t \\
y(t) &= y_0 + v_{0y} t \\
z(t) &= z_0 + v_{0z} t - \tfrac{1}{2} g t^2
\end{aligned}
\qquad
\vec v(t) = [v_{0x},\ v_{0y},\ v_{0z} - g t]
$$

W kodzie robią to funkcje `position` i `velocity`.

### Moment uderzenia

Ciało uderza w powierzchnię, gdy jego wysokość jest równa wysokości powierzchni pod nim:

$$
z(t) = f\big(x(t), y(t)\big)
\quad\Longleftrightarrow\quad
g(t) = z(t) - f\big(x(t), y(t)\big) = 0
$$

Szukamy **najmniejszego dodatniego** $t$, które spełnia to równanie. Ten czas to $t_I$.

### Odbicie

Do odbicia potrzebny jest wektor normalny do powierzchni w punkcie uderzenia.
Dla powierzchni $z = f(x, y)$ wynosi on

$$
\vec N = \left[-\frac{\partial f}{\partial x},\ -\frac{\partial f}{\partial y},\ 1\right],
\qquad
\vec n = \frac{\vec N}{\lVert \vec N \rVert}
$$

Dla paraboloidy $\partial f/\partial x = 2x$ i $\partial f/\partial y = 2y$ (funkcja `grad_f`).

Odbicie działa jak odbicie światła w lustrze. Składowa prędkości wzdłuż normalnej zmienia znak,
a składowa styczna do powierzchni zostaje bez zmian:

$$
\vec v_{\text{po}} = \sqrt{k}\,\Big[\vec v_{\text{przed}} - 2\,(\vec v_{\text{przed}} \cdot \vec n)\,\vec n\Big]
$$

- Dla $k = 1$ długość wektora prędkości się nie zmienia, więc energia kinetyczna też nie (odbicie idealnie sprężyste).
- Dla $k < 1$ prędkość mnożymy przez $\sqrt{k}$. Energia kinetyczna zależy od $v^2$,
  więc po odbiciu wynosi $k \cdot E_\text{kin}$, czyli przy $k = 0.8$ ciało traci 20% energii kinetycznej.

Funkcje `unit_normal` i `reflect` liczą dokładnie te dwa wzory.

## 3. Algorytm krok po kroku

Całość wykonuje funkcja `simulate`, która powtarza poniższą pętlę 5 razy:

```
p, v = punkt startowy, prędkość początkowa
powtórz 5 razy:
    1. znajdź czas uderzenia tI                      (impact_time)
    2. zapisz trajektorię lotu od 0 do tI (do wykresu)
    3. policz punkt uderzenia  pI = position(p, v, tI)
    4. policz prędkość tuż przed uderzeniem  v_przed = velocity(v, tI)
    5. policz energie w punkcie uderzenia i zapisz wiersz tabeli
    6. odbij prędkość względem normalnej  v = reflect(v_przed, n, k)
    7. punkt uderzenia staje się nowym punktem startowym  p = pI
```

Po każdym odbiciu zaczynamy więc nowy rzut ukośny z nowym punktem startowym i nową prędkością.

### Jak szukany jest czas uderzenia

Funkcja $g(t) = z(t) - f(x(t), y(t))$ mówi, czy ciało jest **nad** powierzchnią ($g > 0$), czy **pod** nią ($g < 0$).
W chwili uderzenia $g$ zmienia znak. Program szuka tej zmiany w dwóch krokach:

1. **Przeszukanie zgrubne** (`impact_time`). Oblicza $g(t)$ co 1 ms, od chwili tuż po starcie aż do 100 s,
   i znajduje pierwszy przedział, w którym znak $g$ się zmienia. Wiadomo wtedy, że uderzenie nastąpiło gdzieś w tym przedziale o długości 1 ms.
2. **Bisekcja** (`bisection`). Dzieli ten przedział na pół, sprawdza, w której połówce nadal jest zmiana znaku,
   i powtarza to dalej. Po około 37 podziałach przedział ma szerokość poniżej $10^{-14}$ s,
   więc czas uderzenia jest znany praktycznie dokładnie.

**Dlaczego wyszukiwanie zaczyna się od `T_EPS = 1e-9`, a nie od zera?**
Zaraz po odbiciu ciało leży dokładnie na powierzchni, więc $g(0) = 0$.
Gdyby program sprawdzał też $t = 0$, uznałby, że uderzenie nastąpiło od razu, i utknąłby w miejscu.
Dlatego chwila $t = 0$ jest pomijana. W treści zadania odpowiada to odrzuceniu zerowego rozwiązania.

Jeśli w ciągu 100 s nie ma żadnej zmiany znaku, ciało już nie uderzy w powierzchnię i symulacja się kończy.

### Energie

W każdym punkcie uderzenia ($m = 1$ kg):

$$
E_\text{kin} = \tfrac{1}{2} m v_\text{przed}^2, \qquad
E_\text{pot} = m g z_I, \qquad
E_\text{total} = E_\text{kin} + E_\text{pot}
$$

Energia kinetyczna jest liczona z prędkości **tuż przed** uderzeniem, tak jak w przykładzie z instrukcji.

## 4. Wyniki

Wyniki są takie same jak w przykładzie z instrukcji (do 4 miejsc po przecinku).

### Bez strat ($k = 1$)

| nr | $x_I$ | $y_I$ | $z_I$ | $\Delta t$ | $E_\text{kin}$ | $E_\text{pot}$ | $E_\text{total}$ |
|---|---|---|---|---|---|---|---|
| 1 | 1.1516 | 12.3033 | 152.6967 | 3.0758 | 483.0327 | 1526.9673 | 2010.0000 |
| 2 | 3.7196 | -5.1916 | 40.7882 | 2.5871 | 1602.1179 | 407.8821 | 2010.0000 |
| 3 | -2.6942 | 0.6630 | 7.6981 | 0.5761 | 1933.0189 | 76.9811 | 2010.0000 |
| 4 | 1.1959 | 0.5497 | 1.7324 | 0.1137 | 1992.6760 | 17.3240 | 2010.0000 |
| 5 | -0.5235 | -0.8590 | 1.0119 | 0.0370 | 1999.8808 | 10.1192 | 2010.0000 |


### Ze stratami ($k = 0.8$)

| nr | $x_I$ | $y_I$ | $z_I$ | $\Delta t$ | $E_\text{kin}$ | $E_\text{pot}$ | $E_\text{total}$ |
|---|---|---|---|---|---|---|---|
| 1 | 1.1516 | 12.3033 | 152.6967 | 3.0758 | 483.0327 | 1526.9673 | 2010.0000 |
| 2 | 3.6590 | -4.7787 | 36.2244 | 2.8242 | 1551.1498 | 362.2437 | 1913.3935 |
| 3 | -2.3182 | 0.3989 | 5.5330 | 0.6011 | 1547.8332 | 55.3303 | 1603.1635 |
| 4 | 0.9060 | 0.5192 | 1.0904 | 0.1094 | 1282.6929 | 10.9040 | 1293.5969 |
| 5 | -0.4207 | -0.8360 | 0.8759 | 0.0421 | 1028.2993 | 8.7590 | 1037.0583 |


## 5. Czy energia całkowita jest stała?

**Bez strat: tak.** W każdym punkcie uderzenia $E_\text{total} = 2010$ J, czyli tyle, ile na starcie:
$\tfrac12 \cdot 1 \cdot (2^2 + 4^2) + 1 \cdot 10 \cdot 200 = 10 + 2000 = 2010$ J.
Wynika to z dwóch faktów:

- w czasie lotu działa tylko grawitacja, która jest siłą zachowawczą, więc energia kinetyczna zamienia się w potencjalną i odwrotnie, a suma się nie zmienia,
- odbicie tylko „obraca” wektor prędkości i nie zmienia jego długości, więc nie zmienia też energii kinetycznej.

Widać też, że ciało zsuwa się coraz niżej w miskę: $z_I$ maleje, $E_\text{pot}$ zamienia się w $E_\text{kin}$,
a odstępy między odbiciami $\Delta t$ są coraz krótsze.

**Ze stratami: nie.** Każde odbicie zabiera 20% energii kinetycznej, którą ciało miało tuż przed nim. Na przykład:

$$
2010 - 0.2 \cdot 483.0327 = 2010 - 96.6065 = 1913.3935\ \text{J}
$$

To dokładnie $E_\text{total}$ w drugim wierszu tabeli. Pierwszy wiersz jest taki sam w obu przypadkach,
bo przed pierwszym odbiciem nie było jeszcze żadnych strat.

## 6. Jak uruchomić i co powstaje

Potrzebne są tylko `numpy` i `matplotlib`.

```bash
python lab2.py
```

Program:

- wypisuje obie tabele w konsoli,
- zapisuje je do `results_no_losses.csv` i `results_with_losses.csv`,
- zapisuje wykresy `trajectories_no_losses.png` i `trajectories_with_losses.png` i wyświetla je w oknie.

Na wykresach szara siatka to powierzchnia, kolorowe linie to kolejne trajektorie (1–5),
czerwone kropki to punkty uderzenia, a czarna kropka to punkt startowy.

## 7. Jak zmienić wariant zadania

Wszystkie dane są na początku pliku `lab2.py`:

- `P0`, `V0`: punkt startowy i prędkość początkowa,
- `K_LOSS`: współczynnik strat,
- `f(x, y)`: wzór powierzchni,
- `grad_f(x, y)`: jej pochodne cząstkowe $[\partial f/\partial x,\ \partial f/\partial y]$. Trzeba je policzyć ręcznie i wpisać razem z nowym `f`.

Wzór powierzchni w tytule wykresu (w funkcji `main`) jest wpisany na sztywno, więc przy zmianie `f` trzeba go poprawić również tam.
