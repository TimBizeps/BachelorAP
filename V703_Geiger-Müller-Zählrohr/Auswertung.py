import matplotlib as mpl
mpl.use('pgf')
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
mpl.rcParams.update({
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'pgf.texsystem': 'lualatex',
    'pgf.preamble': r'\usepackage{unicode-math}\usepackage{siunitx}'
})

x, y, z = np.genfromtxt('Daten.txt', unpack=True)
j, k, l = np.genfromtxt('Daten2.txt', unpack=True)
ys = y/60
ks = k/60

print('Messwerte und Fehler:', ys, '±', (np.sqrt(y)/60))

def f(x, a, b):
    return a*x+b

params, covariance = curve_fit(f, j, ks)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

x_plot = np.linspace(300, 710)

np.savetxt("Parameter.txt", np.column_stack([params, errors]))

plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.errorbar(x, ys, yerr=(np.sqrt(y)/60), fmt='ro', label="Messwerte")
plt.xlabel(r'Spannung am Zählrohr in $\si{\volt}$')
plt.ylabel(r'Zählrate in $\si{\per\second}$')
plt.xlim(300, 710)
plt.ylim(420, 560)
plt.legend(loc="best")
plt.title('Charakteristik des Zählrohrs')
plt.tight_layout()
plt.savefig("Plot.pdf")
