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

x, y, z = np.genfromtxt('daten.txt', unpack=True)

R1 = ufloat(48.1, 0.1)
L = ufloat(0.01011, 0.00003)
C = ufloat(0.000000002098, 0.000000000006)

T = 1/(x*1000)
phi = (y/T)*360
print('Phase = ', phi)

plt.plot(x, phi, 'r+', label="Messdaten")
plt.xlabel(r'Frequenz $\nu$ in $\si{\kilo\hertz}$')
plt.ylabel(r'Phasenverschiebung $\phi$ in rad')
plt.yticks([0, 45, 90, 135, 180],
        [r"$0$", r"$\frac{1}{4} \pi$", r"$\frac{1}{2} \pi$", r"$\frac{3}{4} \pi$", r"$\pi$"])
#plt.xscale('log')
#plt.yscale('log')
plt.ylim(0, 180)
plt.xlim(5, 65)
plt.axvline(x=30.8, ymin=0, ymax=0.25, color='b', ls='--', label=r"$\nu_1$ und $\nu_2$", linewidth=1)
plt.axhline(y=45, xmin=0, xmax=0.43, color='b', ls='--', linewidth=1)
plt.axvline(x=38.8, ymin=0, ymax=0.75, color='b', ls='--', linewidth=1)
plt.axvline(x=34.1, ymin=0, ymax=0.50, color='g', ls='--', label=r"$\nu_0$", linewidth=1)
plt.axhline(y=90, xmin=0, xmax=0.485, color='g', ls='--', linewidth=1)
plt.axhline(y=135, xmin=0, xmax=0.563, color='b', ls='--', linewidth=1)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot5.pdf")
