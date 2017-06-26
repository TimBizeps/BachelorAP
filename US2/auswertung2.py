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

t9, u9, hf9, tgc9 = np.genfromtxt('Loch9.txt', unpack=True)

plt.plot(t9, u9, 'b-', label='Amplitude', linewidth=1)
plt.plot(t9, hf9, 'r-', label='HF-Signal', linewidth=1)
plt.xlabel(r'$t$ in $\si{\micro\second}$')
plt.ylabel(r'$U$ in $\si{\volt}$')
plt.xlim(0, 100)
plt.ylim(-1.1, 1.1)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot2.pdf")
