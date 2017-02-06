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
L = 0.023954
C = 0.0000000007932
Csp = 0.000000000028

Ck, A, F = np.genfromtxt('daten.txt', unpack=True)

Ck = unp.uarray(Ck, 0.03*Ck)
print('Ck =', noms(Ck), '±', stds(Ck))
Ck = Ck*10**(-9)
Anzahl = unp.uarray(A, F)
print('Anzahl =', noms(Anzahl), '±', stds(Anzahl))
nup = 1/(2*np.pi*np.sqrt(L*(C+Csp)))
print('nu+ =', nup)
num = 1/(2*np.pi*unp.sqrt(L*((1/((1/C)+(2/Ck)))+Csp)))
print('nu- =', noms(num), '±', stds(num))
