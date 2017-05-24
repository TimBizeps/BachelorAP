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
ladungene = np.array([noms(ladung[6])])
radiuse = np.array([radius[6]])
ladungmea = np.mean(ladungene)
ladungstm = np.std(ladungene, ddof=1) / np.sqrt(len(ladungene))

qkorr = ladung*(1+(B/(p*radius)))**(-(3/2))

qkorrek = np.array([noms(qkorr[6])])
ladungenef = np.array([stds(ladung[6])])
qkorrekf = np.array([stds(qkorr[6])])

print("Radius", radius)
print("Ladung", ladung)
print("Vielfache", n)
print("Elementarladung:", ladungmea, '+/-', ladungstm)
print("Korrigierte Ladung:", qkorr) #durch Korrektur noch kleiner
np.savetxt("radien.txt", np.column_stack([radius]))
np.savetxt("ladungen.txt", np.column_stack([noms(ladung), stds(ladung)]))
np.savetxt("qkorr.txt", np.column_stack([noms(qkorr), stds(qkorr)]))

x = np.linspace(1e-19,3.0e-19,16000)
y = np.zeros(16000)
plt.errorbar(ladungene, radiuse, xerr=ladungenef, fmt='ro', label='Elementarladung nicht korrigiert')
plt.xlabel(r'Ladung $q$ / $\SI{10e-19}{\coulomb}$')
plt.ylabel(r'Tröpfchenradius $r$ / $\si{\micro\meter}$')
plt.legend(loc="best")
plt.title('Nicht korrigierte Ladung')
plt.tight_layout()
plt.savefig("plot_messwerte.pdf")
