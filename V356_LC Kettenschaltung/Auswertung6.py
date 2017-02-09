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

n, U1, U2, U3 = np.genfromtxt('Daten4.txt', unpack=True)

plt.plot(n, U3, 'g-', label=r"$f = \SI{5}{\kilo\hertz}$ (endlose Kette)")
plt.legend(loc="best")
plt.xlabel(r'Kondensator Nr.')
plt.ylabel(r'Spannung in $\si{\volt}$')
plt.tight_layout()
plt.savefig("Plot7.pdf")
