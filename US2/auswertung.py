import matplotlib as mpl
mpl.use('pgf')
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'lualatex',
    'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}'
})

t11, t21, n, t14, t24, sh = np.genfromtxt('Messwerte.txt', unpack=True)
t6, u6, hf6, tgc6 = np.genfromtxt('Loch6.txt', unpack=True)

t41 = 14.9 #micros
t42 = 13.7 #micros

Cluft = 343 #m/s
Cwasser = 1450 #m/s
Cacryl = 2730 #m/s

a =  4.03  #cm
b =  8.04  #cm
c = 15.025 #cm

#1.8 micros korrektur durch geometrische Bestimmung der Dicke des Schutzbereiches
#(Aussendesignal) aus der ersten B-Scan Messung 4.75 cm entspricht 17 micros
#Schicht 0.5 cm Dick => 1.8 micros

#A-Scan
s11 = (t11-1.8)*0.000001*Cacryl/2 #m
s21 = (t21-1.8)*0.000001*Cacryl/2 #m
grl1 = 0.0804-s11-s21 #m

#Auflösevermögen
s41 = (t41-1.8)*0.000001*Cacryl/2 #m
s42 = (t42-1.8)*0.000001*Cacryl/2 #m

#B-Scan 4MHz
s14 = (t14-1.8)*0.000001*Cacryl/2 #m
s24 = (t24-1.8)*0.000001*Cacryl/2 #m
grl4 = 0.0804-s14-s24 #m

#Herz
sl = 7.8 #cm
#nehmen ESV = 0 an
h = sl-sh #cm
th = 17/4.75 * h #micros
Htat = th*0.000001*1450*100 #cm
Hm = np.mean(Htat)
Hf = np.std(Htat, ddof=1) / np.sqrt(len(Htat))
Hzsm = ufloat(Hm, Hf)
print('s11 =', s11)
print('s21 =', s21)
print('Größe1 =', grl1)
print('s41 =', s41)
print('s42 =', s42)
print('s14 =', s14)
print('s24 =', s24)
print('Größe4 =', grl4)
print('h =', h)
print('th =', th)
print('Htat =', Htat)
print('Hm =', Hm, '+-', Hf)
V = 0.5*np.pi*((4.93/2)**2)*Hzsm #cm³
print('EDV =', V)
Vherz = V * 2 #cm³/s
print('Vherz =', Vherz)

plt.plot(t6, u6, 'b-', label='Amplitude', linewidth=1)
plt.plot(t6, hf6, 'r-', label='HF-Signal', linewidth=1)
plt.xlabel(r'$t$ in $\si{\micro\second}$')
plt.ylabel(r'$U$ in $\si{\volt}$')
plt.xlim(0, 100)
plt.ylim(-1.1, 1.1)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot.pdf")
