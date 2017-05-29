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

U3, I3, U4, I4, U5, I5 = np.genfromtxt('messwerte3.txt', unpack=True)
U5log = np.log(U5)
I5log = np.log(I5)

U5fit = np.array([U5log[0], U5log[1], U5log[2], U5log[3], U5log[4], U5log[5], U5log[6], U5log[7], U5log[8], U5log[9], U5log[10], U5log[11], U5log[12], U5log[13], U5log[14], U5log[15], U5log[16], U5log[17], U5log[18], U5log[19], U5log[20], U5log[21], U5log[22]])
I5fit = np.array([I5log[0], I5log[1], I5log[2], I5log[3], I5log[4], I5log[5], I5log[6], I5log[7], I5log[8], I5log[9], I5log[10], I5log[11], I5log[12], I5log[13], I5log[14], I5log[15], I5log[16], I5log[17], I5log[18], I5log[19], I5log[20], I5log[21], I5log[22]])

def f(x, a, b):
    return a*x+b

params, covariance = curve_fit(f, U5fit, I5fit)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

np.savetxt("Parameter.txt", np.column_stack([params, errors]))

x_plot = np.linspace(0, 6)

plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit', linewidth=1)
plt.plot(U5log, I5log, 'gx', label='Kennlinie 5', linewidth=1)
plt.xlabel(r'$U$')
plt.ylabel(r'$I$')
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot2.pdf")
