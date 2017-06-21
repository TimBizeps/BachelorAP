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

f, U = np.genfromtxt('messwerte.txt', unpack=True)
NUsp, Nr3r41, NUbr1, Nr3r42, NUbr2 = np.genfromtxt('messwerte2.txt', unpack=True)
GUsp, Gr3r41, GUbr1, Gr3r42, GUbr2 = np.genfromtxt('messwerte3.txt', unpack=True)
DUsp, Dr3r41, DUbr1, Dr3r42, DUbr2 = np.genfromtxt('messwerte4.txt', unpack=True)

rhowN = 7.24 #g/cm³
rhowG = 7.4 #g/cm³
rhowD = 7.8 #g/cm³

e = const.e #C
me = const.m_e #kg
hq = const.hbar #Js
mü0 = const.mu_0 #T²m³/J
k = const.k #J/K
NA = const.N_A #1/mol

MNd2O3 = 336.48 #g/mol Quelle
MGd2O3 = 362.5 #g/mol  Quelle
MDy2O3 = 373 #g/mol    Quelle
müb = (1/2) * (e/me) * hq #CJs/kg=J/T            1
NN = 1000000 * 2 * NA * rhowN / MNd2O3 #1/m³     2
NG = 1000000 * 2 * NA * rhowG / MGd2O3 #1/m³     2
ND = 1000000 * 2 * NA * rhowD / MDy2O3 #1/m³     2
Ngj = 0.72                                      #4
Ggj = 2                                         #4
Dgj = 1.33                                      #4
Nj = 4.5                                        #3
Gj = 3.5                                        #3
Dj = 7.5                                        #3

Nchitheo = (mü0*(müb**2)*(Ngj**2)*NN*Nj*(Nj+1))/(3*k*293) #293K ~ 20°C
Gchitheo = (mü0*(müb**2)*(Ggj**2)*NG*Gj*(Gj+1))/(3*k*293)
Dchitheo = (mü0*(müb**2)*(Dgj**2)*ND*Dj*(Dj+1))/(3*k*293)

print('müb =', müb)
print('NN =', NN)
print('NG =', NG)
print('ND =', ND)
print('Nchitheo =', Nchitheo)
print('Gchitheo =', Gchitheo)
print('Dchitheo =', Dchitheo)

Nm = 9 #g
Nl = 16.5 #cm
Gm = 14.08 #g
Gl = 16.7 #cm
Dm = 15.1 #g
Dl = 15.8 #cm

n = 250
F = 0.866 #cm²
l = 13.5 #cm
R = 0.7 #Ohm

QrealN = Nm/(Nl*rhowN) #cm²
QrealG = Gm/(Gl*rhowG) #cm²
QrealD = Dm/(Dl*rhowD) #cm²

print('QrealN in cm² =', QrealN)
print('QrealG in cm² =', QrealG)
print('QrealD in cm² =', QrealD)

NUbr = (NUbr2-NUbr1) #mV
GUbr = (GUbr2-GUbr1) #mV
DUbr = (DUbr2-DUbr1) #mV

NUsp = NUsp*1000 #mV
GUsp = GUsp*1000 #mV
DUsp = DUsp*1000 #mV

Nchi1 = 4*(F*NUbr)/(QrealN*NUsp)
Gchi1 = 4*(F*GUbr)/(QrealG*GUsp)
Dchi1 = 4*(F*DUbr)/(QrealD*DUsp)

NdeltaR = Nr3r41 - Nr3r42 #Ohm
GdeltaR = Gr3r41 - Gr3r42 #Ohm
DdeltaR = Dr3r41 - Dr3r42 #Ohm

Nr3 = 998 #Ohm
Gr3 = 998 #Ohm
Dr3 = 998 #Ohm

Nchi2 = 2*(NdeltaR*F)/(Nr3*QrealN)
Gchi2 = 2*(GdeltaR*F)/(Gr3*QrealG)
Dchi2 = 2*(DdeltaR*F)/(Dr3*QrealD)

print('Nchi1 =', Nchi1)
print('Gchi1 =', Gchi1)
print('Dchi1 =', Dchi1)

print('Nchi2 =', Nchi2)
print('Gchi2 =', Gchi2)
print('Dchi2 =', Dchi2)

Nchim1 = np.mean(Nchi1)
Nchif1 = np.std(Nchi1, ddof=1) / np.sqrt(len(Nchi1))
Gchim1 = np.mean(Gchi1)
Gchif1 = np.std(Gchi1, ddof=1) / np.sqrt(len(Gchi1))
Dchim1 = np.mean(Dchi1)
Dchif1 = np.std(Dchi1, ddof=1) / np.sqrt(len(Dchi1))
Nchim2 = np.mean(Nchi2)
Nchif2 = np.std(Nchi2, ddof=1) / np.sqrt(len(Nchi2))
Gchim2 = np.mean(Gchi2)
Gchif2 = np.std(Gchi2, ddof=1) / np.sqrt(len(Gchi2))
Dchim2 = np.mean(Dchi2)
Dchif2 = np.std(Dchi2, ddof=1) / np.sqrt(len(Dchi2))

print('Nchim1 =', Nchim1, '±', Nchif1)
print('Nchim2 =', Nchim2, '±', Nchif2)
print('Gchim1 =', Gchim1, '±', Gchif1)
print('Gchim2 =', Gchim2, '±', Gchif2)
print('Dchim1 =', Dchim1, '±', Dchif1)
print('Dchim2 =', Dchim2, '±', Dchif2)

plt.plot(f, U, 'rx', label='Messwerte', linewidth=1)
plt.xlabel(r'$f$ in $\si{\hertz}$')
plt.ylabel(r'$U_A$ in $\si{\milli\volt}$')
plt.xlim(29900, 40000)
plt.ylim(0, 1050)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot.pdf")
