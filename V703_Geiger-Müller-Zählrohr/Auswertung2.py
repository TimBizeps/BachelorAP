import matplotlib as mpl
mpl.use('pgf')
import numpy as np
from scipy.constants import e
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

x, y, z = np.genfromtxt('Daten.txt', unpack=True)
zs = z/1000000
yw = y/60
deltaQ = zs/yw
elem = deltaQ/e
#def f(x, a, b):
#    return a*x+b
#
#params, covariance = curve_fit(f, x, deltaQ)
#
#errors = np.sqrt(np.diag(covariance))
#
#print('a =', params[0], '±', errors[0])
#print('b =', params[1], '±', errors[1])
#
#x_plot = np.linspace(300, 710)
#
np.savetxt("Ladungen.txt", np.column_stack([deltaQ, elem]))
#
#plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.plot(x, elem, 'rx', label='Messwerte')
plt.xlabel(r'Spannung am Zählrohr in $\si{\volt}$')
plt.ylabel(r'Ladung pro Teilchen in Vielfachen der Elementarladung')
plt.xlim(300, 710)
#plt.ylim(420, 560)
plt.legend(loc="best")
plt.title('Die im Zählrohr freigesetzte Ladung pro Teilchen')
plt.tight_layout()
plt.savefig("Plot2.pdf")
print('e =', e)
