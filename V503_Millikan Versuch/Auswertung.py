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

t_0, U, T, b = np.genfromtxt('V503.txt', unpack=True)
d = ufloat(0.0076250,0.0000051)
s = 0.0005      # [m]
B = 8.22e-3    # [Pa/m]
v = s/t_0
g = const.g     # [m/s^2]
r = 886 #Öldichte
e = const.e
b = b*10**(-5)
p = 101325

radius = np.sqrt((9*b*v)/(2*g*r))
ladung = (4*np.pi/3)*(radius**3)*r*g*d/U
n = ladung/e
ladungene = np.array([noms(ladung[1]), noms(ladung[2]), noms(ladung[4]), noms(ladung[5]), noms(ladung[6]), noms(ladung[13]), noms(ladung[16]), noms(ladung[18]), noms(ladung[19])])
radiuse = np.array([radius[1], radius[2], radius[4], radius[5], radius[6], radius[13], radius[16], radius[18], radius[19]])

qkorr = ladung*(1+(B/(p*radius)))**(-(3/2))

qkorrek = np.array([noms(qkorr[1]), noms(qkorr[2]), noms(qkorr[4]), noms(qkorr[5]), noms(qkorr[6]), noms(qkorr[13]), noms(qkorr[16]),  noms(qkorr[18]), noms(qkorr[19])])

n2 = qkorrek/e

ladungenef = np.array([stds(ladung[1]), stds(ladung[2]), stds(ladung[4]), stds(ladung[5]), stds(ladung[6]), stds(ladung[13]), stds(ladung[16]), stds(ladung[18]), stds(ladung[19])])
qkorrekf = np.array([stds(qkorr[1]), stds(qkorr[2]), stds(qkorr[4]), stds(qkorr[5]), stds(qkorr[6]), stds(qkorr[13]), stds(qkorr[16]), stds(qkorr[18]), stds(qkorr[19])])

faktor = np.array([6, 10, 1, 2, 1, 1, 2, 2, 3])
faktor2 = np.array([5, 8, 1, 2, 1, 1, 1, 2, 2])

elem = ladungene/faktor
elemf = ladungenef/faktor
elemkorr = qkorrek/faktor2
elemkorrf = qkorrekf/faktor2

print("Fehler Elem", elemf)

ladungmea = np.mean(elem)
ladungstm = np.std(elem, ddof=1) / np.sqrt(len(elem))
korrladungmea = np.mean(elemkorr)
korrladungstm = np.std(elemkorr,ddof=1) / np.sqrt(len(elemkorr))

print("Radius", radius)
print("Ladung", ladung)
print("Korrigierte Ladung:", qkorr)
print("elemlad", elem, '+/-', elemf)
print("korrelemlad", elemkorr, '+/-', elemkorrf)
print("Vielfache", n)
print("Vielfache 2", n2)
print("Elementarladung:", ladungmea, '+/-', ladungstm)
print("Korrigierte Elementarladung:", korrladungmea, '+/-', korrladungstm)
np.savetxt("radien.txt", np.column_stack([radius]))
np.savetxt("ladungen.txt", np.column_stack([noms(ladung), stds(ladung)]))
np.savetxt("qkorr.txt", np.column_stack([noms(qkorr), stds(qkorr)]))

x = np.linspace(1e-19,3.0e-19,16000)
y = np.zeros(16000)
plt.errorbar(elemkorr, radiuse, xerr=elemkorrf, fmt='ro', label='Korrigierte Elementarladung')
plt.xlabel(r'Ladung $q$ / $\si{\coulomb}$')
plt.ylabel(r'Tröpfchenradius $r$ / $\si{\micro\meter}$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("plot_messwerte+.pdf")
