# Header
import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit
# from scipy.stats import stats
# import scipy.constants as const
# import scipy.integrate as integrate
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 18


def linear(x, m, b):
    return m*x + b


# Daten einlesen
dt, N = np.genfromtxt('kali.txt', unpack=True)
daten = np.genfromtxt('Daten/Daten.Spe', unpack=True)
kal1 = np.genfromtxt("Daten/Kal_1_micro_second.Spe", unpack=True)
print("Summe Stopp: ", ufloat(np.sum(daten), np.sum(np.sqrt(daten))))
dt = dt     # Zentrieren
N *= 1
hoehe = np.mean(N[4:16])
links = N[0:4]
rechts = N[16:]
dt_rechts = dt[16:]

# linearer Fit links und rechts
params_links, cov_links = curve_fit(linear, dt[0:4], links)
errors_links = np.sqrt(np.diag(cov_links))
m = ufloat(params_links[0], errors_links[0])
b = ufloat(params_links[1], errors_links[1])
print("Steigung links: ", m)
print("y-Achsenabschnitt rechts: ", b)
print("Halbe Höhe: ", hoehe/2)

params_rechts, cov_rechts = curve_fit(linear, dt_rechts, rechts)
errors_rechts = np.sqrt(np.diag(cov_rechts))
m = ufloat(params_rechts[0], errors_rechts[0])
b = ufloat(params_rechts[1], errors_rechts[1])
print("Steigung rechts: ", m)
print("y-Achsenabschnitt rechts: ", b)

# Berechnen des Schnittpunktes
x_links = np.linspace(-25, -17)
x_rechts = np.linspace(6, 20)

links_w = linear(x_links, *params_links)
rechts_w = linear(x_rechts, *params_rechts)

plt.figure(1)
plt.ylabel(r"$N(t) \, / \, (20\mathrm{s})^{-1}$")
plt.xlabel(r"$\mathrm{d}t \, / \, \mathrm{ns}$")
plt.errorbar(dt, N, yerr=np.sqrt(N), fmt='kx', label="Messwerte")
# plt.plot(dt, N, 'r.', label="Messwerte")
plt.axhline(y=hoehe, xmin=0.15, xmax=0.75, label="Plateau")
plt.axhline(y=hoehe/2, xmin=0.065, xmax=0.9, color="green",
            label="Halbwertsbreite")
# plt.ylim(0, 130)
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig("Plateau20.pdf")
plt.clf()

#  Kanalauswertung
kal_t, kanal, hits = np.genfromtxt("Daten/Kalibrierung_1.txt", unpack=True)

params_kal, cov_kal = curve_fit(linear, kanal, kal_t)
errors_kal = np.sqrt(np.diag(cov_kal))
m = ufloat(params_kal[0], errors_kal[0])
b = ufloat(params_kal[1], errors_kal[1])
print("Steigung: ", m)
print("y-Achsenabschnitt: ", b)


# Plot dazu
x = np.linspace(0, 210)
plt.plot(kanal, kal_t, 'r+', label="Daten", markersize=25)
plt.plot(x, linear(x, *params_kal), 'b--', label="Regression")
plt.xlabel("Kanal")
plt.ylabel(r"$T_{VZ} \, / \, \mathrm{\mu s}$")
plt.xlim(0, 210)
plt.tight_layout()
plt.legend(loc="best")
plt.savefig("kal20.pdf")
plt.clf()

# Bestimmung Untergrundrate
messdauer = 147182  # Sekunden
Startimpulse = 3061879
Startimpulse = ufloat(Startimpulse, np.sqrt(Startimpulse))
n = Startimpulse/messdauer
Ts = 20*10**(-6)    # Sekunden
Nf = Startimpulse*n*Ts*unp.exp(-n*Ts)
Nf_kanal = Nf/450
print("-------------------")
print(Startimpulse)
print("Fehlmessungen: ", Nf)
print("Untergrundrate: ", Nf_kanal)

# Umrechnung Kanäle in Zeit
kanaele = np.arange(0, 510, 1)
zeiten1 = linear(kanaele, *params_kal)
# Rausnehmen von komischen Werten
zeiten = zeiten1[0:6]
zeiten = np.append(zeiten, zeiten1[7:14])
zeiten = np.append(zeiten, zeiten1[15:476])

daten_ang = daten[0:6]
daten_ang = np.append(daten_ang, daten[7:14])
daten_ang = np.append(daten_ang, daten[15:4])

# print("Vorher: ", len(daten_ang))
# print("Vorher ohne Nullen: ", len(daten_ang[daten_ang > 0]))

# Nullenverarbeitung
for i in range(len(daten_ang)):
    if daten_ang[i] == 0:
        a = np.array([daten_ang[i-1], daten_ang[i+1]])
        daten_ang[i] = np.round(np.mean(a))

# print("Nachher: ", len(daten_ang))
# print("Nachher ohne Nullen: ", len(daten_ang[daten_ang > 0]))
# print("Länge Zeiten: ", len(zeiten))
# Entfernte Daten und Zeiten
zeiten_ent = zeiten1[:3]
zeiten_ent = np.append(zeiten_ent, zeiten1[6])
zeiten_ent = np.append(zeiten_ent, zeiten1[14])
daten_ent = daten[:3]
daten_ent = np.append(daten_ent, daten[6])
daten_ent = np.append(daten_ent, daten[14])
print(zeiten.shape)
print(daten_ang.shape)
# Definition der exp-Funktion


def e(x, N_0, l, U):
    return N_0*np.exp(-l*x) + U


params_fit, cov_fit = curve_fit(e, zeiten1, daten)
errors_fit = np.sqrt(np.diag(cov_fit))
N_0 = ufloat(params_fit[0], errors_fit[0])
l = ufloat(params_fit[1], errors_fit[1])
U = ufloat(params_fit[2], errors_fit[2])

print("-------------------")
print("N_0: ", N_0)
print("Lambda: ", l)
print("Lebensdauer: ", 1/l)
print("Untergrundrate: ", U)
tau = 2.2
ta_fit = 1/l
print("Verhältnis: ", (1-ta_fit/tau)*100)
print("Verhältnis Untergrund: ", (1-U/Nf_kanal)*100)

print(zeiten1)
plt.hist(zeiten1, bins=100, label="Daten")
# plt.plot(zeiten1, daten, 'r.', label="Daten")
#plt.plot(zeiten_ent, daten_ent, 'gx', label="Nicht-betrachtete Daten")
plt.plot(zeiten, e(zeiten, *params_fit), 'b--', label="Fit")
plt.xlabel(r"$\tau \,  / \, \mathrm{\mu s}$")
plt.ylabel(r"$N(t)$")
plt.grid()
plt.legend(loc="best")
plt.tight_layout()
plt.savefig('fig20.pdf')
plt.clf()


# Tabellen
#np.savetxt('Tabellen/tvz.txt', np.column_stack([dt, N]),
#           delimiter=' & ', newline=r' \\'+'\n', fmt="%.0f")
np.savetxt('Tabellen/kal.txt', np.column_stack([kal_t, kanal, hits]),
           delimiter=' & ', newline=r' \\'+'\n', fmt="%.0f")
np.savetxt('Tabellen/fit.txt', np.column_stack([kanaele, zeiten1, daten]),
           delimiter=' & ', newline=r' \\'+'\n', fmt="%.0f")
print("Alles ohne Probleme ausgeführt!")
