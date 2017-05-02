import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import uncertainties as unp
# Allgemein gilt: N = counts/t

#----------------------------- Aufgabe a) -----------------------------#
#      Bestimmung der Absorptionskoeffizienten und der Größe N(0)      #
#    für Zink und Eisen bei Verwendung eines Cs-137 Gamma-Strahlers    #
#----------------------------------------------------------------------#

Nulleffekt = 1000/1000 # in [Anzahl/Sekunde]

#-------------------------------- Zink --------------------------------#

# Daten laden
d_Pb, t_Pb, counts_Pb = np.loadtxt('V704a.txt', unpack = True)

# In SI-Basiseinheiten umrechnen
d_Pb *= 1e-3

# Auswertung
ln_Pb = np.log(counts_Pb/t_Pb-Nulleffekt)



def f(x, m, b):
    return m*x+b

params_Pb, cov_Pb = curve_fit(f, d_Pb, ln_Pb)
m_Pb  = params_Pb[0]
b_Pb  = params_Pb[1]
Δm_Pb = np.sqrt(cov_Pb[0][0])
Δb_Pb = np.sqrt(cov_Pb[1][1])

print('Geradensteigung Zink:')
print('{}+-{}'.format(m_Pb,Δm_Pb))
print('y-Achsenabschnit Zink:')
print('{}+-{}'.format(b_Pb,Δb_Pb))

x = np.linspace(0, 0.0025, 400)
#plt.plot(d_Pb*1e3, (counts_Pb/t_Pb-Nulleffekt), 'rx', label='Messwerte')
plt.errorbar(d_Pb*1e3, (counts_Pb/t_Pb-Nulleffekt), np.sqrt(counts_Pb/t_Pb-Nulleffekt), fmt='rx', label='Messwerte')
plt.plot(x*1e3, np.exp(f(x, m_Pb, b_Pb)), 'k-', label='Ausgleichsfunktion')
plt.yscale('log')
#plt.xlim(0, 2.5)
#plt.ylim(10, 250)
#plt.ylim(2, 5.5)
plt.xlabel(r'$\mathrm{Absorberschichtdicke}\; d/\mathrm{mm}$')
plt.ylabel(r'$\mathrm{Strahlintensität}\; N-N_{\mathrm{u}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_Zink.pdf")
plt.close()

#-------------------------------- Eisen -------------------------------#

# Daten laden
d_Fe, t_Fe, counts_Fe = np.loadtxt('V704b.txt', unpack = True)

# In SI-Basiseinheiten umrechnen
d_Fe *= 1e-3

# Auswertung
ln_Fe = np.log(counts_Fe/t_Fe-Nulleffekt)

def f(x, m, b):
    return m*x+b

params_Fe, cov_Fe = curve_fit(f, d_Fe, ln_Fe)
m_Fe  = params_Fe[0]
b_Fe  = params_Fe[1]
Δm_Fe = np.sqrt(cov_Fe[0][0])
Δb_Fe = np.sqrt(cov_Fe[1][1])

print(counts_Fe/t_Fe-Nulleffekt)
print(ln_Fe)

print('Geradensteigung Eisen:')
print('{}+-{}'.format(m_Fe,Δm_Fe))
print('y-Achsenabschnit Eisen:')
print('{}+-{}'.format(b_Fe,Δb_Fe))

x = np.linspace(0, 0.005, 400)
#plt.plot(d_Fe*1e3, (counts_Fe/t_Fe-Nulleffekt), 'rx', label='Messwerte')
plt.errorbar(d_Fe*1e3, (counts_Fe/t_Fe-Nulleffekt), np.sqrt(counts_Fe/t_Fe-Nulleffekt), fmt='rx', label='Messwerte')
plt.plot(x*1e3, np.exp(f(x, m_Fe, b_Fe)), 'k-', label='Ausgleichsfunktion')
plt.yscale('log')
#plt.xlim(-0.2, 2.5)
#plt.ylim(10, 250)
#plt.ylim(2, 5.5)
plt.xlabel(r'$\mathrm{Absorberschichtdicke}\; d/\mathrm{mm}$')
plt.ylabel(r'$\mathrm{Strahlintensität}\; N-N_{\mathrm{u}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_Eisen.pdf")
plt.close()

#----------------------------- Aufgabe b) -----------------------------#
#                                                                      #
#----------------------------------------------------------------------#

Nulleffekt = 327/1000 # in [Anzahl/Sekunde]

d_A, d_B, t_A, N_A = np.loadtxt('V704c.txt', unpack = True)

d_A *= 1e-6

ln_A = np.log(N_A/t_A-Nulleffekt)

def f(x, m, b):
    return m*x+b

params_A1, cov_A1 = curve_fit(f, d_A[0:5], ln_A[0:5])
m_A1  = params_A1[0]
b_A1  = params_A1[1]
Δm_A1 = np.sqrt(cov_A1[0][0])
Δb_A1 = np.sqrt(cov_A1[1][1])

def g(x, b):
    return b

params_A2, cov_A2 = curve_fit(g, d_A[5:10], ln_A[5:10])
#m_A2  = params_A2[0]
b_A2  = params_A2[0]
#Δm_A2 = np.sqrt(cov_A2[0][0])
Δb_A2 = np.sqrt(cov_A2[0][0])




m_A2 = 0
#b_A2 = np.log(0.1876967839)
Δm_A2 = 0
#Δb_A2 = 0

print('Geradensteigung 1 Alu:')
print('{}+-{}'.format(m_A1,Δm_A1))
print('y-Achsenabschnit 1 Alu:')
print('{}+-{}'.format(b_A1,Δb_A1))

print('Geradensteigung 2 Alu:')
print('{}+-{}'.format(m_A2,Δm_A2))
print('y-Achsenabschnit 2 Alu:')
print('{}+-{}'.format(b_A2,Δb_A2))

R = (b_A2-b_A1)/(m_A1-m_A2)
ΔR = np.sqrt((1/(m_A1-m_A2)*Δb_A2)**2+(1/(m_A1-m_A2)*Δb_A1)**2+((b_A1-b_A2)/(m_A2-m_A1)**2*Δm_A1)**2+((b_A2-b_A1)/(m_A1-m_A2**2)*Δm_A2)**2)

rho=2.7 # g/cm^3
E = 1.92*np.sqrt((rho*R*1e2)**2+0.22*(rho*R*1e2))
ΔE = 1.92*np.sqrt((rho*ΔR*1e2)**2+0.22*(rho*ΔR*1e2))

print(R)
print(ΔR)
print(E)
print(ΔE)

x = np.linspace(100e-6, 500e-6, 2000)
#plt.plot(d_A*1e6, N_A/t_A-Nulleffekt, 'rx', label='Messwerte')
plt.errorbar(d_A*1e6, N_A/t_A-Nulleffekt, np.sqrt(N_A/t_A-Nulleffekt), fmt='rx', label='Messwerte')
plt.plot(x*1e6, np.exp(f(x, m_A1, b_A1)), 'k-', label='Ausgleichsfunktionen')
plt.plot(x*1e6, np.exp(f(x, m_A2, b_A2)), 'k-')
plt.yscale('log', nonposy='clip')
#plt.xlim(100, 500)
#plt.ylim(0.1, 100)
plt.xlabel(r'$\mathrm{Absorberschichtdicke}\; d/\mathrm{\mu m}$')
plt.ylabel(r'$\mathrm{Strahlintensität}\; N-N_{\mathrm{u}}$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig("plot_Aluminium.pdf")
plt.close()
