import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
from scipy.stats import stats
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18

# Erster Teil
# Einlesen der Daten
pumpPro15, fmeanHz15, inten15 = np.genfromtxt('winkel15.txt', unpack=True)
pumpPro30, fmeanHz30, inten30 = np.genfromtxt('winkel30.txt', unpack=True)
pumpPro60, fmeanHz60, inten60 = np.genfromtxt('winkel60.txt', unpack=True)

# Bestimmung der Momentangeschwindigkeit
f0Hz = 2 * 10 ** 6
c = 1800
v15 = fmeanHz15*c/(2*f0Hz*10) * 1/np.cos(80.1)
v30 = fmeanHz30*c/(2*f0Hz) * 1/np.cos(70.5)
v60 = fmeanHz60*c/(2*f0Hz) * 1/np.cos(54.7)
print(1/np.cos(80.1))
print(pumpPro60)
print(fmeanHz15)
print("15 iz da:")
print(v15)
print(pumpPro60)
print(fmeanHz30)
print("30 iz da:")
print(v30)
print(pumpPro60)
print(fmeanHz60)
print("60 iz da:")
print(v60)

# Plots
plt.figure(1)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v15, fmeanHz15/np.cos(80.1) * 10 ** (-3), 'ro', label="Dopplerwinkel 80.1°")
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a15.pdf')

plt.figure(2)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v30, fmeanHz30/np.cos(70.5) * 10 ** (-3), 'go', label="Dopplerwinkel 70.5°")
#plt.xticks([0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6],
#            [])
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a30.pdf')

plt.figure(3)
plt.xlabel(r"$v/ \mathrm{\frac{m}{s}}$")
plt.ylabel(r"$\frac{\mathrm{\Delta} \nu}{\mathrm{cos} \alpha} / \mathrm{kHz}$")
plt.plot(v60, fmeanHz60/np.cos(54.7) * 10 ** (-3), 'bo', label="Dopplerwinkel 54.7°")
plt.legend(loc="best")
plt.tight_layout
plt.savefig('a60.pdf')

# Aufgabenteil 2
# Auslesen der Daten
tiefeMikroS45, fmeanHz45, inten45 = np.genfromtxt('teilb45.txt', unpack=True)
tiefeMikroS70, fmeanHz70, inten70 = np.genfromtxt('teilb70.txt', unpack=True)

# Tiefe umrechnen
tiefeMM45 = tiefeMikroS45 * 1.5 + 17.25
tiefeMM70 = tiefeMikroS70 * 1.5 + 17.25

# Bestimmung der Momentangeschwindigkeit
v45 = fmeanHz45*c/(2*f0Hz) * 1/np.cos(80.1)
v70 = fmeanHz70*c/(2*f0Hz) * 1/np.cos(80.1)
print(tiefeMM45)
print("45er:")
print(v45)
print(inten45)
print(tiefeMM70)
print("70er:")
print(v70)
print(inten70)

# Plots
plt.figure(4)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel(r"$\mathrm{Messtiefe} / \mathrm{mm}$")
ax1.set_ylabel(r"$I / \mathrm{\frac{V^2}{s} * 10}$")
ax2.set_ylabel(r"$v / \mathrm{\frac{m}{s}}$")
ax1.plot(tiefeMM45, inten45 *10 ** (-2), 'r+', label="$I$")
ax2.plot(tiefeMM45, v45, 'bo', label="$v$")
ax1.set_xlim(35, 48)
ax1.set_ylim(-15, 10)
ax2.set_ylim(-15, 10)
ax1.legend(loc="lower left")
ax2.legend(loc="lower right")
plt.tight_layout
plt.savefig('b45.pdf')

plt.figure(5)
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.set_xlabel(r"$\mathrm{Messtiefe} / \mathrm{mm}$")
ax1.set_ylabel(r"$I / \mathrm{\frac{V^2}{s} * 10}$")
ax2.set_ylabel(r"$v / \mathrm{\frac{m}{s}}$")
ax1.plot(tiefeMM70, inten70 *10 ** (-2), 'r+', label="$I$")
ax2.plot(tiefeMM70, v70, 'bo', label="$v$")
ax1.set_xlim(35, 48)
ax1.set_ylim(-30, 70)
ax2.set_ylim(-30, 70)
ax1.legend(loc="lower left")
ax2.legend(loc="lower right")
plt.tight_layout
plt.savefig('b70.pdf')
