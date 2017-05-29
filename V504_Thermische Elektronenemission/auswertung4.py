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

I1 = 1.8 #A
U1 = 3.25 #V
Is1 = 0.000029 #A
I2 = 1.9 #A
U2 = 3.5 #V
Is2 = 0.000068 #A
I3 = 2 #A
U3 = 4 #V
Is3 = 0.000156 #A
I4 = 2.1 #A
U4 = 4.25 #V
Is4 = 0.000340 #A
I5 = 2.2 #A
U5 = 4.5 #V
Is5 = 0.000715 #A

e = const.e #C
m = const.m_e #kg
h = const.h #Js
k = const.k #J/K
sigma = 5.7*(10**-12) #W/cm²K⁴
eta = 0.28
fcm = 0.32 #cm²
fm = 0.000032 #m²

print('e', e)
print('k', k)
print('h', h)
print('m', m)

Nwl = 0.95 #W

N1 = U1*I1 #W
N2 = U2*I2 #W
N3 = U3*I3 #W
N4 = U4*I4 #W
N5 = U5*I5 #W

print('N1 =', N1)
print('N2 =', N2)
print('N3 =', N3)
print('N4 =', N4)
print('N5 =', N5)

T1 = ((N1 - Nwl)/(fcm*eta*sigma))**(1/4) #K
T2 = ((N2 - Nwl)/(fcm*eta*sigma))**(1/4) #K
T3 = ((N3 - Nwl)/(fcm*eta*sigma))**(1/4) #K
T4 = ((N4 - Nwl)/(fcm*eta*sigma))**(1/4) #K
T5 = ((N5 - Nwl)/(fcm*eta*sigma))**(1/4) #K

print('T1 =', T1)
print('T2 =', T2)
print('T3 =', T3)
print('T4 =', T4)
print('T5 =', T5)

js1 = Is1/fm #A/m²
js2 = Is2/fm #A/m²
js3 = Is3/fm #A/m²
js4 = Is4/fm #A/m²
js5 = Is5/fm #A/m²

efi1 = - np.log((js1*(h**3))/(4*np.pi*e*m*(k**2)*(T1**2)))*k*T1 #Austrittsarbeit J
efi2 = - np.log((js2*(h**3))/(4*np.pi*e*m*(k**2)*(T2**2)))*k*T2 #Austrittsarbeit J
efi3 = - np.log((js3*(h**3))/(4*np.pi*e*m*(k**2)*(T3**2)))*k*T3 #Austrittsarbeit J
efi4 = - np.log((js4*(h**3))/(4*np.pi*e*m*(k**2)*(T4**2)))*k*T4 #Austrittsarbeit J
efi5 = - np.log((js5*(h**3))/(4*np.pi*e*m*(k**2)*(T5**2)))*k*T5 #Austrittsarbeit J

efi1 = efi1*6.242*10**18 #in eV
efi2 = efi2*6.242*10**18 #in eV
efi3 = efi3*6.242*10**18 #in eV
efi4 = efi4*6.242*10**18 #in eV
efi5 = efi5*6.242*10**18 #in eV

print('js1 =', js1)
print('js2 =', js2)
print('js3 =', js3)
print('js4 =', js4)
print('js5 =', js5)

print('efi1 =', efi1)
print('efi2 =', efi2)
print('efi3 =', efi3)
print('efi4 =', efi4)
print('efi5 =', efi5)

efi = np.array([efi1, efi2, efi3, efi4, efi5])

efim = np.mean(efi)
efistm = np.std(efi, ddof=1) / np.sqrt(len(efi))

print('efi =', efim, '±', efistm)
