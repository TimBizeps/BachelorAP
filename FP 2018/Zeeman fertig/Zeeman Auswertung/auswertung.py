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

DeltaSr, deltas = np.genfromtxt('daten3.txt', unpack=True)

quoti = deltas/DeltaSr

DeltaLambDr = 4.89*10**(-11)
DeltaLambDb = 2.7*10**(-11)

print('rot ', quoti)

DeltaSb, deltas1, deltas2, deltas31, deltas32, deltas4 = np.genfromtxt('daten4.txt', unpack=True)

quoti1 = deltas1/DeltaSb
quoti2 = deltas2/DeltaSb
quoti31 = deltas31/DeltaSb
quoti32 = deltas32/DeltaSb
quoti4 = deltas4/DeltaSb

print('b1 ', quoti1)
print('b2 ', quoti2)
print('b31 ', quoti31)
print('b32 ', quoti32)
print('b4 ', quoti4)

DeltaSb2, deltas5 = np.genfromtxt('daten5.txt', unpack=True)

quoti5 = deltas5/DeltaSb2

print('b5 ', quoti5)

deltaLambr = (1/2)*quoti*DeltaLambDr
deltaLambb1 = (1/2)*quoti1*DeltaLambDb
deltaLambb2 = (1/2)*quoti2*DeltaLambDb
deltaLambb31 = (1/2)*quoti31*DeltaLambDb
deltaLambb32 = (1/2)*quoti32*DeltaLambDb
deltaLambb4 = (1/2)*quoti4*DeltaLambDb
deltaLambb5 = (1/2)*quoti5*DeltaLambDb
print('deltaLambr ', deltaLambr)
print('deltaLambb1 ', deltaLambb1)
print('deltaLambb2 ', deltaLambb2)
print('deltaLambb31 ', deltaLambb31)
print('deltaLambb32 ', deltaLambb32)
print('deltaLambb4 ', deltaLambb4)
print('deltaLambb5 ', deltaLambb5)
deltaLambb3 = np.concatenate((deltaLambb31, deltaLambb32), axis=None)
print('deltaLambb3 ', deltaLambb3)
meandeltalar = np.mean(deltaLambr)
fdeltalar = np.std(deltaLambr, ddof=1)/np.sqrt(len(deltaLambr))
meandeltalab1 = np.mean(deltaLambb1)
fdeltalab1 = np.std(deltaLambb1, ddof=1)/np.sqrt(len(deltaLambb1))
meandeltalab2 = np.mean(deltaLambb2)
fdeltalab2 = np.std(deltaLambb2, ddof=1)/np.sqrt(len(deltaLambb2))
meandeltalab3 = np.mean(deltaLambb3)
fdeltalab3 = np.std(deltaLambb3, ddof=1)/np.sqrt(len(deltaLambb3))
meandeltalab4 = np.mean(deltaLambb4)
fdeltalab4 = np.std(deltaLambb4, ddof=1)/np.sqrt(len(deltaLambb4))
meandeltalab5 = np.mean(deltaLambb5)
fdeltalab5 = np.std(deltaLambb5, ddof=1)/np.sqrt(len(deltaLambb5))

deltaLaR = ufloat(meandeltalar, fdeltalar)
deltaLaB1 = ufloat(meandeltalab1, fdeltalab1)
deltaLaB2 = ufloat(meandeltalab2, fdeltalab2)
deltaLaB3 = ufloat(meandeltalab3, fdeltalab3)
deltaLaB4 = ufloat(meandeltalab4, fdeltalab4)
deltaLaB5 = ufloat(meandeltalab5, fdeltalab5)

print('deltaLaR ', deltaLaR)
print('deltaLaB1 ', deltaLaB1)
print('deltaLaB2 ', deltaLaB2)
print('deltaLaB3 ', deltaLaB3)
print('deltaLaB4 ', deltaLaB4)
print('deltaLaB5 ', deltaLaB5)

#h =
#c =
#müB = 

#E = (h*c)/Lambda
#E = gj*müB*B
#=> gj = (h*c)/(Lambda*müB*B)
