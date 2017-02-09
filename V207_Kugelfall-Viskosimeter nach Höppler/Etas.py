import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp

#Dichten der Kugeln
d_kl = ufloat(1.564, 0.005)
m_kl = ufloat(4.46, 0.01)
dichte_kl = m_kl/((4/3)*np.pi*((d_kl/2)**3))

d_gr = ufloat(1.582, 0.005)
m_gr = ufloat(4.60, 0.01)
dichte_gr = m_gr/((4/3)*np.pi*((d_gr/2)**3))

np.savetxt("Dichten.txt", np.column_stack([dichte_kl.n, dichte_kl.s, dichte_gr.n, dichte_gr.s]))

dichte_H2O_20 = 0.9982
dichte_H2O_25 = 0.9971
dichte_H2O_30 = 0.9957
dichte_H2O_35 = 0.9940
dichte_H2O_40 = 0.9932
dichte_H2O_45 = 0.9902
dichte_H2O_50 = 0.9880
dichte_H2O_55 = 0.9857
dichte_H2O_60 = 0.9832
dichte_H2O_65 = 0.9806
dichte_H2O_70 = 0.9778

K_kl = 0.0007640

#Zeiten
a, b = np.genfromtxt('DatenRaumT.txt', unpack=True)

t_gr = np.mean(a)
delta_t_gr = 1/np.sqrt(len(a)) * np.std(a, ddof=1)
t_kl = np.mean(b)
delta_t_kl = 1/np.sqrt(len(b)) * np.std(b, ddof=1)

c, d, e, f, g, h, i, j, k, l = np.genfromtxt('DatenTemperaturen.txt', unpack=True)

t_gr_25 = np.mean(c)
delta_t_25 = 1/np.sqrt(len(c)) * np.std(c, ddof=1)
t_gr_30 = np.mean(d)
delta_t_30 = 1/np.sqrt(len(d)) * np.std(d, ddof=1)
t_gr_35 = np.mean(e)
delta_t_35 = 1/np.sqrt(len(e)) * np.std(e, ddof=1)
t_gr_40 = np.mean(f)
delta_t_40 = 1/np.sqrt(len(f)) * np.std(f, ddof=1)
t_gr_45 = np.mean(g)
delta_t_45 = 1/np.sqrt(len(g)) * np.std(g, ddof=1)
t_gr_50 = np.mean(h)
delta_t_50 = 1/np.sqrt(len(h)) * np.std(h, ddof=1)
t_gr_55 = np.mean(i)
delta_t_55 = 1/np.sqrt(len(i)) * np.std(i, ddof=1)
t_gr_60 = np.mean(j)
delta_t_60 = 1/np.sqrt(len(j)) * np.std(j, ddof=1)
t_gr_65 = np.mean(k)
delta_t_65 = 1/np.sqrt(len(k)) * np.std(k, ddof=1)
t_gr_70 = np.mean(l)
delta_t_70 = 1/np.sqrt(len(l)) * np.std(l, ddof=1)

#Zeiten als UFloats
t_k = ufloat(t_kl, delta_t_kl)
t_g = ufloat(t_gr, delta_t_gr)
t_g25 = ufloat(t_gr_25, delta_t_25)
t_g30 = ufloat(t_gr_30, delta_t_30)
t_g35 = ufloat(t_gr_35, delta_t_35)
t_g40 = ufloat(t_gr_40, delta_t_40)
t_g45 = ufloat(t_gr_45, delta_t_45)
t_g50 = ufloat(t_gr_50, delta_t_50)
t_g55 = ufloat(t_gr_55, delta_t_55)
t_g60 = ufloat(t_gr_60, delta_t_60)
t_g65 = ufloat(t_gr_65, delta_t_65)
t_g70 = ufloat(t_gr_70, delta_t_70)

ZeitenN = [t_kl, t_gr, t_gr_25, t_gr_30, t_gr_35, t_gr_40, t_gr_45, t_gr_50, t_gr_55, t_gr_60, t_gr_65, t_gr_70]
ZeitenS = [delta_t_kl, delta_t_gr, delta_t_25, delta_t_30, delta_t_35, delta_t_40, delta_t_45, delta_t_50, delta_t_55, delta_t_60, delta_t_65, delta_t_70]

np.savetxt("Zeiten.txt", np.column_stack([ZeitenN, ZeitenS]))

eta_20 = K_kl * (dichte_kl-dichte_H2O_20) * t_k

#K der großen Kugel
K_gr = eta_20/((dichte_gr-dichte_H2O_20) * t_g)

print(K_gr)

#etas(T)
eta_25 = K_gr * (dichte_gr-dichte_H2O_25) * t_g25
eta_30 = K_gr * (dichte_gr-dichte_H2O_30) * t_g30
eta_35 = K_gr * (dichte_gr-dichte_H2O_35) * t_g35
eta_40 = K_gr * (dichte_gr-dichte_H2O_40) * t_g40
eta_45 = K_gr * (dichte_gr-dichte_H2O_45) * t_g45
eta_50 = K_gr * (dichte_gr-dichte_H2O_50) * t_g50
eta_55 = K_gr * (dichte_gr-dichte_H2O_55) * t_g55
eta_60 = K_gr * (dichte_gr-dichte_H2O_60) * t_g60
eta_65 = K_gr * (dichte_gr-dichte_H2O_65) * t_g65
eta_70 = K_gr * (dichte_gr-dichte_H2O_70) * t_g70

lneta_20 = unp.log(eta_20)
lneta_25 = unp.log(eta_25)
lneta_30 = unp.log(eta_30)
lneta_35 = unp.log(eta_35)
lneta_40 = unp.log(eta_40)
lneta_45 = unp.log(eta_45)
lneta_50 = unp.log(eta_50)
lneta_55 = unp.log(eta_55)
lneta_60 = unp.log(eta_60)
lneta_65 = unp.log(eta_65)
lneta_70 = unp.log(eta_70)

EtasN = [unp.nominal_values(lneta_20), unp.nominal_values(lneta_25), unp.nominal_values(lneta_30), unp.nominal_values(lneta_35), unp.nominal_values(lneta_40), unp.nominal_values(lneta_45), unp.nominal_values(lneta_50), unp.nominal_values(lneta_55), unp.nominal_values(lneta_60), unp.nominal_values(lneta_65), unp.nominal_values(lneta_70)]
EtasS = [unp.std_devs(lneta_20), unp.std_devs(lneta_25), unp.std_devs(lneta_30), unp.std_devs(lneta_35), unp.std_devs(lneta_40), unp.std_devs(lneta_45), unp.std_devs(lneta_50), unp.std_devs(lneta_55), unp.std_devs(lneta_60), unp.std_devs(lneta_65), unp.std_devs(lneta_70)]

np.savetxt("DatenPlot.txt", np.column_stack([EtasN, EtasS]))

v_20_kl=10/t_k
v_20=10/t_g
v_25=10/t_g25
v_30=10/t_g30
v_35=10/t_g35
v_40=10/t_g40
v_45=10/t_g45
v_50=10/t_g50
v_55=10/t_g55
v_60=10/t_g60
v_65=10/t_g65
v_70=10/t_g70

L = d_gr
Re_20_kl = (v_20_kl*d_kl*dichte_H2O_20)/eta_20
Re_20 = (v_20*L*dichte_H2O_20)/eta_20
Re_25 = (v_25*L*dichte_H2O_25)/eta_25
Re_30 = (v_30*L*dichte_H2O_30)/eta_30
Re_35 = (v_35*L*dichte_H2O_35)/eta_35
Re_40 = (v_40*L*dichte_H2O_40)/eta_40
Re_45 = (v_45*L*dichte_H2O_45)/eta_45
Re_50 = (v_50*L*dichte_H2O_50)/eta_50
Re_55 = (v_55*L*dichte_H2O_55)/eta_55
Re_60 = (v_60*L*dichte_H2O_60)/eta_60
Re_65 = (v_65*L*dichte_H2O_65)/eta_65
Re_70 = (v_70*L*dichte_H2O_70)/eta_70

ReynoldsZ = [Re_20_kl.n, Re_20.n, Re_25.n, Re_30.n, Re_35.n, Re_40.n, Re_45.n, Re_50.n, Re_55.n, Re_60.n, Re_65.n, Re_70.n]
ReynoldsZFehler = [Re_20_kl.s, Re_20.s, Re_25.s, Re_30.s, Re_35.s, Re_40.s, Re_45.s, Re_50.s, Re_55.s, Re_60.s, Re_65.s, Re_70.s]

np.savetxt("ReynoldscheZahlen.txt", np.column_stack([ReynoldsZ, ReynoldsZFehler]))
