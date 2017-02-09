import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

x, y = np.genfromtxt('daten.txt', unpack=True)

def f(x, m, b):
    return m * x + b

params, covariance = curve_fit(f, x, y)
# covariance is the covariance matrix

errors = np.sqrt(np.diag(covariance))

print('m =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])

x_plot = np.linspace(0, 0.045)

plt.plot(x, y, 'rx', label="example data")
plt.plot(x_plot, f(x_plot, *params), 'b-', label='linearer Fit', linewidth=3)
plt.legend(loc="best")
plt.xlabel(r'$T^2 \,/\, s^2$')
plt.ylabel(r'$a^2 \,/\, m^2$')
plt.savefig("Plot.pdf")

m_zm=0.44589
m=ufloat(params[0], errors[0])
b=ufloat(params[1], errors[1])

D=(4*np.pi**2*m_zm)/m

print('D =', D)

I_0 = 0.00298775
I_d = (b*D/(4*np.pi**2))-I_0

print('I_0 =', I_0)
print('I_d =', I_d)

T_Zylinder = ufloat(1.19, 0.0001666)
T_Kugel = ufloat(1.697, 0.003679)
T_PuppeI = ufloat(0.8575, 0.0035)
T_PuppeII = ufloat(0.65625, 0.00263)

I_Zylinder = ((T_Zylinder**2*D)/(4*np.pi**2))-I_d
I_Kugel = ((T_Kugel**2*D)/(4*np.pi**2))-I_d
I_PuppeI = ((T_PuppeI**2*D)/(4*np.pi**2))-I_d
I_PuppeII = ((T_PuppeII**2*D)/(4*np.pi**2))-I_d

print('I_Zylinder =', I_Zylinder)
print('I_Kugel =', I_Kugel)
print('I_PuppeI =', I_PuppeI)
print('I_PuppeII =', I_PuppeII)
