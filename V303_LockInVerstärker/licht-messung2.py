import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit

d, U = np.genfromtxt('messwerte2.txt', unpack=True)

d=d-4.9
d=d/100

logd = np.log(d)
logU = np.log(U)

plt.plot(d, U, 'k.', label='$Abstandsverhältnis der Lichtintensität$')
plt.xlabel(r'$\mathrm{d}/m$')
plt.xscale('log')
plt.ylabel(r'$\mathrm{U}/V$')
plt.yscale('log')
#plt.ylim(0.02, 10)

def f(x, b, a):
    return

params, cov = curve_fit(f, d, U)

errors= np.sqrt(np.diag(cov))
print('			Parameter für Data3 gefittet:	')
print('			b=', params[0], '\pm', errors[0])
print('			a=', params[1], '\pm', errors[1])

x_plot = np.linspace(0.05, 1.50)
# plt.plot(d, U, 'r.', label= 'Messwerte')
plt.plot(x_plot, f(x_plot, *params), 'b-', label='Fit')

plt.tight_layout()
plt.savefig('licht-messung2.pdf')
# plt.show()
