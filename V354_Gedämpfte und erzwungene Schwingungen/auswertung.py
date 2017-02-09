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

#x, y, z = np.genfromtxt('daten.txt', unpack=True)

R1 = ufloat(48.1, 0.1)
R2 = ufloat(509.5, 0.5)
L = ufloat(0.01011, 0.00003)
C = ufloat(0.000000002098, 0.000000000006)

q = 1/R2 * unp.sqrt(L/C)
print('q = ', noms(q),'', stds(q))
tex = (2*L)/(R1)
print('tex = ', noms(tex),'', stds(tex))
nüres = (unp.sqrt((1/(L*C))-((R2**2)/(2*L**2))))/(2*np.pi)
print('nüres = ', noms(nüres),'', stds(nüres))
nünü = R2/(L*2*np.pi)
print('nünü = ', noms(nünü),'', stds(nünü))
nü = unp.sqrt(1/(L*C))/(2*np.pi)
print('nü0 = ', noms(nü),'', stds(nü))
na = (-(R2/(2*L)) + unp.sqrt(((R2**2)/(4*(L**2)))+(1/(L*C))))/(2*np.pi)
print('nü1 = ', noms(na),'', stds(na))
no = ((R2/(2*L)) + unp.sqrt(((R2**2)/(4*(L**2)))+(1/(L*C))))/(2*np.pi)
print('nü2 = ', noms(no),'', stds(no))
print(np.pi)
#plt.plot(x, z, 'r+', label="Messdaten")
#plt.xlabel(r'Frequenz $\nu$ in $\si{\kilo\hertz}$')
#plt.ylabel(r'Amplitude $\frac{U_C}{U_0}$')
##plt.xscale('log')
##plt.yscale('log')
#plt.xlim(0, 50)
#plt.ylim(0, 4)
#plt.axvline(x=28, ymin=0, ymax=0.654, color='b', ls='--', label=r"$\nu_+$ und $\nu_-$", linewidth=1)
#plt.axhline(y=3.7, xmin=0, xmax=0.67, color='g', ls='--', linewidth=1)
#plt.axvline(x=38, ymin=0, ymax=0.654, color='b', ls='--', linewidth=1)
#plt.axvline(x=33.5, ymin=0, ymax=0.925, color='g', ls='--', label="Resonanzmaximum", linewidth=1)
#plt.axhline(y=2.616, xmin=0.58, xmax=0.76, color='r', ls='--', label="Resonanzbreite", linewidth=1)
#plt.legend(loc="best")
#plt.tight_layout()
#plt.savefig("Plot3.pdf")
