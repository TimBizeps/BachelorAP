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
#erst gerechnet, dann gemittelt
r21, r31, r41 = np.genfromtxt('daten.txt', unpack=True)
lin1 = r31/r41
rx1 = r21*lin1
Rx1m = np.mean(rx1)
Rx1f = np.std(rx1, ddof=1)/np.sqrt(len(rx1))
print('Rx1m = ', Rx1m)
print('Rx1f = ', Rx1f)
r21f = 0.002*r21
lin1f = 0.005*lin1
R21 = unp.uarray(r21, r21f)
Lin1 = unp.uarray(lin1, lin1f)
Rx1 = R21*Lin1
print('Rx1 = ', noms(Rx1), '±', stds(Rx1))
stdsrx1 = np.array(stds(Rx1))
Rx1f = np.sqrt((stdsrx1[0]/3)**2+(stdsrx1[1]/3)**2+(stdsrx1[2]/3)**2)
print('Rx1f = ', Rx1f)

r22, r32, r42 = np.genfromtxt('daten2.txt', unpack=True)
lin2 = r32/r42
rx2 = r22*lin2
Rx2m = np.mean(rx2)
Rx2f = np.std(rx2, ddof=1)/np.sqrt(len(rx2))
print('Rx2m = ', Rx2m)
print('Rx2f = ', Rx2f)
r22f = 0.002*r22
lin2f = 0.005*lin2
R22 = unp.uarray(r22, r22f)
Lin2 = unp.uarray(lin2, lin2f)
Rx2 = R22*Lin2
print('Rx2 = ', noms(Rx2), '±', stds(Rx2))
stdsrx2 = np.array(stds(Rx2))
Rx2f = np.sqrt((stdsrx2[0]/3)**2+(stdsrx2[1]/3)**2+(stdsrx2[2]/3)**2)
print('Rx2f = ', Rx2f)

c23, r33, r43 = np.genfromtxt('daten3.txt', unpack=True)
lin3 = r33/r43
cx1 = c23*(1/lin3)
Cx1m = np.mean(cx1)
Cx1f = np.std(cx1, ddof=1)/np.sqrt(len(cx1))
print('Cx1m = ', Cx1m)
print('Cx1f = ', Cx1f)
c23f = 0.002*c23
lin3f = 0.005*lin3
C23 = unp.uarray(c23, c23f)
Lin3 = unp.uarray(lin3, lin3f)
Cx1 = C23*(1/Lin3)
print('Cx1 = ', noms(Cx1), '±', stds(Cx1))
stdscx1 = np.array(stds(Cx1))
Cx1f = np.sqrt((stdscx1[0]/3)**2+(stdscx1[1]/3)**2+(stdscx1[2]/3)**2)
print('Cx1f = ', Cx1f)

c24, r34, r44 = np.genfromtxt('daten4.txt', unpack=True)
lin4 = r34/r44
cx2 = c24*(1/lin4)
Cx2m = np.mean(cx2)
Cx2f = np.std(cx2, ddof=1)/np.sqrt(len(cx2))
print('Cx2m = ', Cx2m)
print('Cx2f = ', Cx2f)
c24f = 0.002*c24
lin4f = 0.005*lin4
C24 = unp.uarray(c24, c24f)
Lin4 = unp.uarray(lin4, lin4f)
Cx2 = C24*(1/Lin4)
print('Cx2 = ', noms(Cx2), '±', stds(Cx2))
stdscx2 = np.array(stds(Cx2))
Cx2f = np.sqrt((stdscx2[0]/3)**2+(stdscx2[1]/3)**2+(stdscx2[2]/3)**2)
print('Cx2f = ', Cx2f)

c25, r25, r35, r45, l26, r26, r36, r46, r27, r37, r47 = np.genfromtxt('daten5.txt', unpack=True)
lin5 = r35/r45
cx3 = c25*(1/lin5)
rx3 = r25*lin5
Cx3m = np.mean(cx3)
Cx3f = np.std(cx3, ddof=1)/np.sqrt(len(cx3))
print('Cx3m = ', Cx3m)
print('Cx3f = ', Cx3f)
Rx3m = np.mean(rx3)
Rx3f = np.std(rx3, ddof=1)/np.sqrt(len(rx3))
print('Rx3m = ', Rx3m)
print('Rx3f = ', Rx3f)
c25f = 0.002*c25
r25f = 0.03*r25
lin5f = 0.005*lin5
C25 = unp.uarray(c25, c25f)
R25 = unp.uarray(r25, r25f)
Lin5 = unp.uarray(lin5, lin5f)
Cx3 = C25*(1/Lin5)
Rx3 = R25*Lin5
print('Cx3 = ', noms(Cx3), '±', stds(Cx3))
print('Rx3 = ', noms(Rx3), '±', stds(Rx3))
stdscx3 = np.array(stds(Cx3))
Cx3f = np.sqrt((stdscx3[0]/3)**2+(stdscx3[1]/3)**2+(stdscx3[2]/3)**2)
print('Cx3f = ', Cx3f)
stdsrx3 = np.array(stds(Rx3))
Rx3f = np.sqrt((stdsrx3[0]/3)**2+(stdsrx3[1]/3)**2+(stdsrx3[2]/3)**2)
print('Rx3f = ', Rx3f)

lin6 = r36/r46
lx1 = l26*lin6
rxl1 = r26*lin6
Lx1m = np.mean(lx1)
Lx1f = np.std(lx1, ddof=1)/np.sqrt(len(lx1))
print('Lx1m = ', Lx1m)
print('Lx1f = ', Lx1f)
Rxl1m = np.mean(rxl1)
Rxl1f = np.std(rxl1, ddof=1)/np.sqrt(len(rxl1))
print('Rxl1m = ', Rxl1m)
print('Rxl1f = ', Rxl1f)
l26f = 0.002*l26
r26f = 0.03*r26
lin6f = 0.005*lin6
L26 = unp.uarray(l26, l26f)
R26 = unp.uarray(r26, r26f)
Lin6 = unp.uarray(lin6, lin6f)
Lx1 = L26*Lin6
Rxl1 = R26*Lin6
print('Lx1 = ', noms(Lx1), '±', stds(Lx1))
print('Rxl1 = ', noms(Rxl1), '±', stds(Rxl1))
stdslx1 = np.array(stds(Lx1))
Lx1f = np.sqrt((stdslx1[0]/3)**2+(stdslx1[1]/3)**2+(stdslx1[2]/3)**2)
print('Lx1f = ', Lx1f)
stdsrxl1 = np.array(stds(Rxl1))
Rxl1f = np.sqrt((stdsrxl1[0]/3)**2+(stdsrxl1[1]/3)**2+(stdsrxl1[2]/3)**2)
print('Rxl1f = ', Rxl1f)

c4 = 0.000000994 #nF
rxl2 = r27*r37/r47
lx2 = r27*r37*c4
Lx2m = np.mean(lx2)
Lx2f = np.std(lx2, ddof=1)/np.sqrt(len(lx2))
print('Lx2m = ', Lx2m)
print('Lx2f = ', Lx2f)
Rxl2m = np.mean(rxl2)
Rxl2f = np.std(rxl2, ddof=1)/np.sqrt(len(rxl2))
print('Rxl2m = ', Rxl2m)
print('Rxl2f = ', Rxl2f)
c4f = 0.002*c4
r27f = 0.002*r27
r37f = 0.03*r37
r47f = 0.03*r47
C4 = unp.uarray(c4, c4f)
R27 = unp.uarray(r27, r27f)
R37 = unp.uarray(r37, r37f)
R47 = unp.uarray(r47, r47f)
Rxl2 = R27*R37/R47
Lx2 = R27*R37*C4
print('Lx2 = ', noms(Lx2), '±', stds(Lx2))
print('Rxl2 = ', noms(Rxl2), '±', stds(Rxl2))
stdslx2 = np.array(stds(Lx2))
Lx2f = np.sqrt((stdslx2[0]/3)**2+(stdslx2[1]/3)**2+(stdslx2[2]/3)**2)
print('Lx2f = ', Lx2f)
stdsrxl2 = np.array(stds(Rxl2))
Rxl2f = np.sqrt((stdsrxl2[0]/3)**2+(stdsrxl2[1]/3)**2+(stdsrxl2[2]/3)**2)
print('Rxl2f = ', Rxl2f)


f, Us, Ubr = np.genfromtxt('daten6.txt', unpack=True)
Ubrampl = Ubr/2
Ubreff = Ubrampl/np.sqrt(2)
omega0 = 1/((Cx2m/1000000000)*1000) #R = 1kOhm und C = Cx (Wert3)
Frequenz = omega0/(2*np.pi)
print('Omega0 = ', omega0)
print('Frequenz = ', Frequenz)
print('Ubreff = ', Ubreff)
nü = f
nü0 = 378

y = Ubreff/Us
Omega = nü/nü0
print('y = ', y)
print('Omega = ', Omega)

def f(x):
    return np.sqrt((1/9)*((x**2-1)**2/(((1-x**2)**2)+9*x**2)))

x_plot=np.linspace(0.7, 50, 1000)

plt.plot(Omega, y, 'rx', label="Messwerte")
plt.plot(x_plot, f(x_plot), 'b-', label='Vergleichsfunktion', linewidth=1)
plt.xscale('log')
plt.legend(loc="best")
plt.xlabel(r'Kreisfrequenzverhältniss: $\frac{\nu}{\nu_0}$')
plt.ylabel(r'Spannungsverhältniss: $\frac{U_\mathup{Br_{eff}}}{U_\mathup{S}}$')
plt.tight_layout()
plt.savefig("Plot.pdf")

Ubrnü0 = Ubr[0]
Ubrnü0eff = Ubrnü0/(2*np.sqrt(2))
U1 = Us[0]
U2 = Ubrnü0eff/(np.sqrt((2**2-1)**2/(9*((1-2**2)**2+9*2**2))))#Formel vom Plot mit Omega 2

k=U2/U1
print('Ubrnü0eff = ', Ubrnü0eff)
print('U2 = ', U2)
print('k = ', k)
