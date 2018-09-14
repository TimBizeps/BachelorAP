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

r1 = ufloat(470, 4)
rp = ufloat(33100, 200)
Ub1 = ufloat(13.0, 0.1)
Ub2 = ufloat(-14.0, 0.1)

Uee = (r1/rp)*(-Ub2)
print('Uee = ', Uee)

Uea = -(r1/rp)*Ub1
print('Uea = ', Uea)
