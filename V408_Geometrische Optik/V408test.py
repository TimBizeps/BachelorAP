import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

V, b = np.loadtxt("V408f.txt", unpack=True) # lädt Messdaten aus txt in Variablen g und b_1

def f(x, m, b):
		return m*x + b      # Regressionsfunktion. Hier: linear
params, cov = curve_fit(f, 1/V, b) # fit der Variablen g und b_1 an Funktion f
m = params[0]       # Nenne ersten Parameter des fits m
c = params[1]       # Nenne zweiten Parameter des fits b
dm = np.sqrt(cov[0][0]) # aus der Kovarianzmatrix lassen sich die Fehler des fits ablesen
db = np.sqrt(cov[1][1]) # dito
print('''
----------------------------------------------------
Steigung der Regression: {:5f}+-{}

Achsenabschnitt der Regression: {:5f}+-{}
----------------------------------------------------'''.format(m, dm, c, db))
print(b)
print(V)
x = np.linspace(0, 4, 10000)  # erzeuge ganz viele x Werte um fit zu plotten
plt.plot(x, f(x,m,c), 'k-', label='Regression')
plt.plot(1/V, b, 'rx', label='Messwerte')
# plt.yscale('log') Y-Achse logarithmisieren
# plt.xlim(-0.2, 2.5) X-Achse Begrenzungen
# plt.ylim(10, 250) Y-Achse Begrenzungen
plt.xlabel(r'$\mathrm{V}$')
plt.ylabel(r'$\mathrm{Bildweite\,\,b^\prime}\; \mathrm{[cm]}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("AbbeV352363.pdf")
plt.close()
