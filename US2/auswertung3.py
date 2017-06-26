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

def f(x):
    return 0.54*x

x_plot = np.linspace(0, 15)

plt.plot(x_plot, f(x_plot), 'b-', linewidth=1)
plt.ylabel(r'$\si{\centi\meter}$')
plt.xlim(0, 15)
plt.ylim(0, 8)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot3.png")
