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

m, f = np.genfromtxt('Daten3.txt', unpack=True)

def v(omega, L, C):
    return omega/np.arccos(1-(1/2)*(omega**2)*L*C)

l = 0.00175
c1 = 0.000000022
c2 = 0.00000000939
phase = m * np.pi/16
kreisfreq = 2 * np.pi * f
phasegesch = kreisfreq/phase

x_plot = np.linspace(100, 245000)
plt.plot(kreisfreq, phasegesch, 'rx', label="Messdaten")
plt.plot(x_plot, v(x_plot, 0.00175, 0.000000022), 'b-', label='Theoriekurve', linewidth=1)
plt.legend(loc="best")
plt.xticks([0, 50000, 100000, 150000, 200000, 250000, 300000, 350000],
        [0, 50, 100, 150, 200, 250, 300, 350])
plt.yticks([110000, 120000, 130000, 140000, 150000, 160000, 170000],
        [110, 120, 130, 140, 150, 160, 170])
plt.xlabel(r'$\omega$ in $\si{\per\milli\per\second}$')
plt.ylabel(r'$v_\mathup{Ph}$ in $\si{\kilo\per\second}$')
plt.tight_layout()
plt.savefig("Plot6.pdf")
print(phasegesch)
