import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
x = np.array([30, 32, 35, 40])
y = np. array([100.08, 107.34, 118.29, 136.32])

def linear(x,a,b):
    return a*x+b

params, cov = curve_fit(linear, x, y)

m = params[0]
m_err = np.sqrt(cov[0][0])

b = params[1]
b_err = np.sqrt(cov[1][1])

c = np.linspace(0, 50)

plt.plot(x,y, 'rx', label='Messwerte')
plt.plot(c, linear(c,m,b), label='Regression')
plt.grid()
plt.xlabel('Ordnungszahl des Elements')
plt.ylabel('E[eV]')
plt.legend(loc='best')
plt.savefig("plot.pdf")
plt.close()

print(b)
print(b_err)
print(m)
print(m_err)
