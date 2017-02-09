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

u, w, y, z = np.genfromtxt('daten.txt', unpack=True)
u = u + 273.15
w = w + 273.15
def T(z, a, b, c, j):
    return a*z**3+b*z**2+c*z+j

def S(z, d, e, f):
    return ((d*z**2)/(1+(e*z**2)))+f

params, covariance = curve_fit(T, z, u)

errors = np.sqrt(np.diag(covariance))

print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
print('j =', params[3], '±', errors[3])

A=ufloat(params[0],errors[0])
B=ufloat(params[1],errors[1])
C=ufloat(params[2],errors[2])
J=ufloat(params[3],errors[3])

x_plot = np.linspace(0, 1200)

plt.plot(z, u, 'r+', label=r"Messdaten $T_1$")
plt.plot(x_plot, T(x_plot, *params), 'r-', label='Ausgleichskurve', linewidth=1)

params2, covariance2 = curve_fit(S, z, w)

errors2 = np.sqrt(np.diag(covariance2))

print('a =', params2[0], '±', errors2[0])
print('b =', params2[1], '±', errors2[1])
print('c =', params2[2], '±', errors2[2])

D=ufloat(params2[0],errors2[0])
E=ufloat(params2[1],errors2[1])
F=ufloat(params2[2],errors2[2])

plt.plot(z, w, 'b+', label=r"Messdaten $T_2$")
plt.plot(x_plot, S(x_plot, *params2), 'b-', label='Ausgleichskurve', linewidth=1)
plt.xlabel(r'Zeit $t$ in $\si{\second}$')
plt.ylabel(r'Temperatur in $\si{\kelvin}$')
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plot.pdf")

ablt11 = 3*A*120**2+2*B*120+C #3.
ablt12 = 3*A*420**2+2*B*420+C #8.
ablt13 = 3*A*720**2+2*B*720+C #13.
ablt14 = 3*A*960**2+2*B*960+C #17.
ablt21 = (2*D*120)/((1+E*120**2)**2) #3.
ablt22 = (2*D*420)/((1+E*420**2)**2) #8.
ablt23 = (2*D*720)/((1+E*720**2)**2) #13.
ablt24 = (2*D*960)/((1+E*960**2)**2) #17.

wert11 = A*120**3+B*120**2+C*120+J
wert12 = A*420**3+B*420**2+C*420+J
wert13 = A*720**3+B*720**2+C*720+J
wert14 = A*960**3+B*960**2+C*960+J
wert21 = ((D*120**2)/(1+(E*120**2)))+F
wert22 = ((D*420**2)/(1+(E*420**2)))+F
wert23 = ((D*720**2)/(1+(E*720**2)))+F
wert24 = ((D*960**2)/(1+(E*960**2)))+F

print('ablt11 =', ablt11.n, '±', ablt11.s)
print('ablt12 =', ablt12.n, '±', ablt12.s)
print('ablt13 =', ablt13.n, '±', ablt13.s)
print('ablt14 =', ablt14.n, '±', ablt14.s)
print('ablt21 =', ablt21.n, '±', ablt21.s)
print('ablt22 =', ablt22.n, '±', ablt22.s)
print('ablt23 =', ablt23.n, '±', ablt23.s)
print('ablt24 =', ablt24.n, '±', ablt24.s)
#print('wert11 =', wert11.n, '±', wert11.s)
#print('wert12 =', wert12.n, '±', wert12.s)
#print('wert13 =', wert13.n, '±', wert13.s)
#print('wert14 =', wert14.n, '±', wert14.s)
#print('wert21 =', wert21.n, '±', wert21.s)
#print('wert22 =', wert22.n, '±', wert22.s)
#print('wert23 =', wert23.n, '±', wert23.s)
#print('wert24 =', wert24.n, '±', wert24.s)

mkck=660#J/K
m1cw=3*4180#kg*J/(kg*K)

ablq11 = (m1cw+mkck)*ablt11
ablq12 = (m1cw+mkck)*ablt12
ablq13 = (m1cw+mkck)*ablt13
ablq14 = (m1cw+mkck)*ablt14
ablq21 = (m1cw+mkck)*ablt21
ablq22 = (m1cw+mkck)*ablt22
ablq23 = (m1cw+mkck)*ablt23
ablq24 = (m1cw+mkck)*ablt24

print('ablq11 =', ablq11.n, '±', ablq11.s)
print('ablq12 =', ablq12.n, '±', ablq12.s)
print('ablq13 =', ablq13.n, '±', ablq13.s)
print('ablq14 =', ablq14.n, '±', ablq14.s)
print('ablq21 =', ablq21.n, '±', ablq21.s)
print('ablq22 =', ablq22.n, '±', ablq22.s)
print('ablq23 =', ablq23.n, '±', ablq23.s)
print('ablq24 =', ablq24.n, '±', ablq24.s)

N1 = 175.0
N2 = 205.0
N3 = 212.5
N4 = 212.5

g1 = ablq11/N1
g2 = ablq12/N2
g3 = ablq13/N3
g4 = ablq14/N4

print('Güte1 =', g1.n, '±', g1.s)
print('Güte2 =', g2.n, '±', g2.s)
print('Güte3 =', g3.n, '±', g3.s)
print('Güte4 =', g4.n, '±', g4.s)

idg1 = 296.15/(296.15-294.15)
idg2 = 304.95/(304.95-286.15)
idg3 = 314.15/(314.15-277.55)
idg4 = 320.35/(320.35-273.15)

print('ideale Güte1 =', idg1)
print('ideale Güte2 =', idg2)
print('ideale Güte3 =', idg3)
print('ideale Güte4 =', idg4)

v, x, n, o = np.genfromtxt('daten2.txt', unpack=True)

n = n + 273.15
o = o + 273.15
p0 = 1.013
v = v/p0
def R(n, h, i):
    return h * n + i

params3, covariance3 = curve_fit(R, 1/n, np.log(v))

errors3 = np.sqrt(np.diag(covariance3))

print('h =', params3[0], '±', errors3[0])
print('i =', params3[1], '±', errors3[1])

R = 8.314 #joule/(mol*K) Quelle!!!
Steigung = ufloat(params3[0], errors3[0])
L = -Steigung*R
print('L =', L.n, '±', L.s)

mdurchs1 = ablq21/L
mdurchs2 = ablq22/L
mdurchs3 = ablq23/L
mdurchs4 = ablq24/L

print('Massendurchsatz1 =', mdurchs1.n, '±', mdurchs1.s)#Einheit mol/s
print('Massendurchsatz2 =', mdurchs2.n, '±', mdurchs2.s)
print('Massendurchsatz3 =', mdurchs3.n, '±', mdurchs3.s)
print('Massendurchsatz4 =', mdurchs4.n, '±', mdurchs4.s)

Mcl2f2c = 120.9 #g/mol

Mgdurchs1 = mdurchs1 * Mcl2f2c
Mgdurchs2 = mdurchs2 * Mcl2f2c
Mgdurchs3 = mdurchs3 * Mcl2f2c
Mgdurchs4 = mdurchs4 * Mcl2f2c

print('Massendurchsatz1 in g/s =', Mgdurchs1.n, '±', Mgdurchs1.s)#Einheit g/s
print('Massendurchsatz2 in g/s =', Mgdurchs2.n, '±', Mgdurchs2.s)
print('Massendurchsatz3 in g/s =', Mgdurchs3.n, '±', Mgdurchs3.s)
print('Massendurchsatz4 in g/s =', Mgdurchs4.n, '±', Mgdurchs4.s)

rho0 = 5510 #g/m³ = 5.51 g/l
T0 = 273.15 #K
ptr0 = 1 #Bar
kapa = 1.14

rho120 = ((2.6*T0)/(ptr0*294.15))*rho0
rho420 = ((3.2*T0)/(ptr0*286.15))*rho0
rho720 = ((3.2*T0)/(ptr0*277.55))*rho0
rho960 = ((3.2*T0)/(ptr0*273.15))*rho0

print('Dichte 120 =', rho120)#Einheit g/m³
print('Dichte 420 =', rho420)
print('Dichte 720 =', rho720)
print('Dichte 960 =', rho960)

Nmech1 = (1/(kapa-1))*(7*(2.6/7)**(1/1.14)-2.6)*100000*(1/rho120)*Mgdurchs1 #*100000 da 1 Bar = 100000 N/m²
Nmech2 = (1/(kapa-1))*(9*(3.2/9)**(1/1.14)-3.2)*100000*(1/rho420)*Mgdurchs2
Nmech3 = (1/(kapa-1))*(11*(3.2/11)**(1/1.14)-3.2)*100000*(1/rho720)*Mgdurchs3
Nmech4 = (1/(kapa-1))*(12.5*(3.2/12.5)**(1/1.14)-3.2)*100000*(1/rho960)*Mgdurchs4

print('Mechanische Leistung =', Nmech1.n, '±', Nmech1.s)#Einheit W
print('Mechanische Leistung =', Nmech2.n, '±', Nmech2.s)
print('Mechanische Leistung =', Nmech3.n, '±', Nmech3.s)
print('Mechanische Leistung =', Nmech4.n, '±', Nmech4.s)

wirk1 = 100*Nmech1/N1
wirk2 = 100*Nmech2/N2
wirk3 = 100*Nmech3/N3
wirk4 = 100*Nmech4/N4

print('Wirkungsgrad =', wirk1.n, '±', wirk1.s)#Einheit W
print('Wirkungsgrad =', wirk2.n, '±', wirk2.s)
print('Wirkungsgrad =', wirk3.n, '±', wirk3.s)
print('Wirkungsgrad =', wirk4.n, '±', wirk4.s)
