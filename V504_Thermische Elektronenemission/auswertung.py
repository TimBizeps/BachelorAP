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

U1, I1 = np.genfromtxt('messwerte.txt', unpack=True)
U2, I2 = np.genfromtxt('messwerte2.txt', unpack=True)
U3, I3, U4, I4, U5, I5 = np.genfromtxt('messwerte3.txt', unpack=True)

#x_plot = np.linspace(, , )

#plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.plot(U1, I1, 'kx', label='Kennlinie 1', linewidth=1)
plt.plot(U2, I2, 'rx', label='Kennlinie 2', linewidth=1)
plt.plot(U3, I3, 'mx', label='Kennlinie 3', linewidth=1)
plt.plot(U4, I4, 'bx', label='Kennlinie 4', linewidth=1)
plt.plot(U5, I5, 'gx', label='Kennlinie 5', linewidth=1)
plt.yticks([0, 0.029, 0.068, 0.1, 0.156, 0.2, 0.3, 0.340, 0.4, 0.5, 0.6, 0.7, 0.715, 0.8],
        [r"$0.0$", r"$I_\mathup{S_1}$", r"$I_\mathup{S_2}$", r"$0.1$", r"$I_\mathup{S_3}$", r"$0.2$", r"$0.3$", r"$I_\mathup{S_4}$", r"$0.4$", r"$0.5$", r"$0.6$", r"$0.7$", r"$I_\mathup{S_5}$", r"$0.8$"])
plt.xlabel(r'$U$ / $\si{\volt}$')
plt.ylabel(r'$I$ / $\si{\milli\ampere}$')
plt.xlim(-5, 255)
plt.axhline(y=0.715, xmin=0, xmax=1, color='g', ls='-', linewidth=1)
plt.axhline(y=0.340, xmin=0, xmax=1, color='b', ls='-', linewidth=1)
plt.axhline(y=0.156, xmin=0, xmax=1, color='m', ls='-', linewidth=1)
plt.axhline(y=0.068, xmin=0, xmax=1, color='r', ls='-', linewidth=1)
plt.axhline(y=0.029, xmin=0, xmax=1, color='k', ls='-', linewidth=1)
plt.legend(loc="best")
plt.title('Kennlinien')
plt.tight_layout()
plt.savefig("Plot.pdf")
