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

f, Ua = np.genfromtxt('daten1.txt', unpack=True)

plt.loglog(f, Ua, 'rx', label='Messwerte')
plt.xlabel(r'Frequenz $f$ in $\si{\kilo\hertz}$')
plt.ylabel(r'Ausgangsspannung $U_A$ in $\si{\milli\volt}$')
#plt.xlim(,)
plt.ylim(100, 400)
plt.axhline(y=183.5, xmin=0, xmax=1, color='b', ls='--', label=r"$\frac{V^\prime}{\sqrt{2}}$", linewidth=1)
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot.pdf")

#f_g etwa 1150 kHz
