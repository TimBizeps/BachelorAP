import matplotlib as mpl
mpl.use('pgf')
import numpy as np
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
m, f = np.genfromtxt('Daten5.txt', unpack=True)

def omega1(theta, L, C1, C2):
    return np.sqrt((1/(L*C1))+(1/(L*C2))-(1/L)*np.sqrt((((1/C1)+(1/C2))**2)-(4*(np.sin(theta)**2))/(C1*C2)))

def omega2(theta, L, C1, C2):
    return np.sqrt((1/(L*C1))+(1/(L*C2))+(1/L)*np.sqrt((((1/C1)+(1/C2))**2)-(4*(np.sin(theta)**2))/(C1*C2)))

print(omega1(0, 0.00175, 0.000000022, 0.00000000939))

l = 0.00175
c1 = 0.000000022
c2 = 0.00000000939
phase = m * np.pi/16
kreisfreq = 2 * np.pi * f

x_plot = np.linspace(0, np.pi/2)
plt.plot(phase, kreisfreq, 'rx', label="Messdaten")
plt.plot(x_plot, omega1(x_plot, 0.00175, 0.000000022, 0.00000000939), 'b-', label='Theoriekurve: akustischer Ast', linewidth=1)
plt.plot(x_plot, omega2(x_plot, 0.00175, 0.000000022, 0.00000000939), 'g-', label='Theoriekurve: optischer Ast', linewidth=1)
plt.legend(loc="best")
plt.yticks([0, 50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000],
        [0, 50, 100, 150, 200, 250, 300, 350, 400, 450])
plt.xlabel(r'$\theta$ pro Kettenglied in rad')
plt.ylabel(r'$\omega$ in $\si{\per\milli\per\second}$')
plt.tight_layout()
plt.savefig("Plot4.pdf")
