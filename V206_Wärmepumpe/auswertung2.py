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

v, x, n, o = np.genfromtxt('daten2.txt', unpack=True)
n = n + 273.15
o = o + 273.15
p0 = 1.013
v = v/p0
def R(n, h, i):
    return h * n + i

params3, covariance3 = curve_fit(R, 1/n, np.log(v))

errors3 = np.sqrt(np.diag(covariance3))

print('a =', params3[0], '±', errors3[0])
print('b =', params3[1], '±', errors3[1])

x_plot = np.linspace(0.00305, 0.00345)
plt.plot(1/n, np.log(v), 'rx', label=r"Messwerte")
plt.plot(x_plot, R(x_plot, *params3), 'b-', label='Ausgleichsgerade', linewidth=1)
plt.legend(loc="best")
plt.xlabel(r'$\frac{1}{T} \,\,/\,\, \frac{1}{K}$')
plt.ylabel(r'$\ln \left( \frac{p_\mathup{b}}{p_0} \right)$')
plt.tight_layout()
plt.savefig("Plot2.pdf")
