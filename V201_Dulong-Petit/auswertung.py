import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import ufloat
import uncertainties.unumpy as unp
import sympy

U_100 = ufloat(4.10, 0.01)
U_0 = ufloat(-0.002, 0.001)

a = (U_100 - U_0)/100 #U=a*T => T = U/a

Ux = ufloat(0.859, 0.003)
Uy = ufloat(4.10, 0.01)
Uxy = ufloat(2.150, 0.005)

up, vp, wp, xp, yp, zp, uc, vc, wc, xc, yc, zc = np.genfromtxt('Messdaten.txt', unpack=True)
Upw1 = ufloat(up[0], vp[0])
Upk1 = ufloat(wp[0], xp[0])
Upm1 = ufloat(yp[0], zp[0])
Upw2 = ufloat(up[1], vp[1])
Upk2 = ufloat(wp[1], xp[1])
Upm2 = ufloat(yp[1], zp[1])
Upw3 = ufloat(up[2], vp[2])
Upk3 = ufloat(wp[2], xp[2])
Upm3 = ufloat(yp[2], zp[2])

Ucw1 = ufloat(uc[0], vc[0])
Uck1 = ufloat(wc[0], xc[0])
Ucm1 = ufloat(yc[0], zc[0])
Ucw2 = ufloat(uc[1], vc[1])
Uck2 = ufloat(wc[1], xc[1])
Ucm2 = ufloat(yc[1], zc[1])
Ucw3 = ufloat(uc[2], vc[2])
Uck3 = ufloat(wc[2], xc[2])
Ucm3 = ufloat(yc[2], zc[2])

Uaw = ufloat(0.868, 0.001)
Uak = ufloat(3.8, 0.1)
Uam = ufloat(0.988, 0.001)

Tx = (Ux-U_0)/a
Ty = (Uy-U_0)/a
Txy = (Uxy-U_0)/a

Tpw1 = (Upw1-U_0)/a
Tpk1 = (Upk1-U_0)/a
Tpm1 = (Upm1-U_0)/a
Tpw2 = (Upw2-U_0)/a
Tpk2 = (Upk2-U_0)/a
Tpm2 = (Upm2-U_0)/a
Tpw3 = (Upw3-U_0)/a
Tpk3 = (Upk3-U_0)/a
Tpm3 = (Upm3-U_0)/a

Tcw1 = (Ucw1-U_0)/a
Tck1 = (Uck1-U_0)/a
Tcm1 = (Ucm1-U_0)/a
Tcw2 = (Ucw2-U_0)/a
Tck2 = (Uck2-U_0)/a
Tcm2 = (Ucm2-U_0)/a
Tcw3 = (Ucw3-U_0)/a
Tck3 = (Uck3-U_0)/a
Tcm3 = (Ucm3-U_0)/a

Taw = (Uaw-U_0)/a
Tak = (Uak-U_0)/a
Tam = (Uam-U_0)/a

print('a =', a)
print('1/a =', 1/a)
print('Tx =', Tx)
print('Ty =', Ty)
print('Txy =', Txy)
print('Tpw1 =', Tpw1)
print('Tpk1 =', Tpk1)
print('Tpm1 =', Tpm1)
print('Tpw2 =', Tpw2)
print('Tpk2 =', Tpk2)
print('Tpm2 =', Tpm2)
print('Tpw3 =', Tpw3)
print('Tpk3 =', Tpk3)
print('Tpm3 =', Tpm3)
print('Tcw1 =', Tcw1)
print('Tck1 =', Tck1)
print('Tcm1 =', Tcm1)
print('Tcw2 =', Tcw2)
print('Tck2 =', Tck2)
print('Tcm2 =', Tcm2)
print('Tcw3 =', Tcw3)
print('Tck3 =', Tck3)
print('Tcm3 =', Tcm3)
print('Taw =', Taw)
print('Tak =', Tak)
print('Tam =', Tam)

mp = 385.76
mc = 114.55
ma = 108.72
mx = 346.27          #Wichtig für die Rechnung cgmg
my_vorher = 300.23   #Berechning von mx
my_nachher = 271.67  #Wichtig für die Rechnung cgmg
mxy_vorher = 646.50  #Berechnung von mx
mxy_nachher = 617.94 #Wichtig für die Rechnung cgmg

#Massen mw
mwp1 = 593.97
mwp2 = 575.55
mwp3 = 566.31
mwc1 = 580.33
mwc2 = 592.24
mwc3 = 578.50
mwa = 580.05

cw = 4.18

cgmg = (cw*my_nachher*(Ty-Txy)-cw*mx*(Txy-Tx))/(Txy-Tx)

print('cgmg =', cgmg)

cp1 = ((cw*mwp1 + cgmg)*(Tpm1-Tpw1))/(mp*(Tpk1-Tpm1))
cp2 = ((cw*mwp2 + cgmg)*(Tpm2-Tpw2))/(mp*(Tpk2-Tpm2))
cp3 = ((cw*mwp3 + cgmg)*(Tpm3-Tpw3))/(mp*(Tpk3-Tpm3))
cc1 = ((cw*mwc1 + cgmg)*(Tcm1-Tcw1))/(mc*(Tck1-Tcm1))
cc2 = ((cw*mwc2 + cgmg)*(Tcm2-Tcw2))/(mc*(Tck2-Tcm2))
cc3 = ((cw*mwc3 + cgmg)*(Tcm3-Tcw3))/(mc*(Tck3-Tcm3))
ca = ((cw*mwa + cgmg)*(Tam-Taw))/(ma*(Tak-Tam))

print('cp1 =', cp1)
print('cp2 =', cp2)
print('cp3 =', cp3)
print('cc1 =', cc1)
print('cc2 =', cc2)
print('cc3 =', cc3)
print('ca =', ca)

#Cp
Mp = 207.2
rhop = 11350000
alphap = 29*10**(-6)
kp = 42000000000
Mc = 12.0
rhoc = 2250000
alphac = 8*10**(-6)
kc = 33000000000
Ma = 27.0
rhoa = 2700000
alphaa = 23.5*10**(-6)
ka = 75000000000
V0p = Mp/rhop
V0c = Mc/rhoc
V0a = Ma/rhoa

R = 8.314

Cpq_p_1 = cp1*Mp
Cpq_p_2 = cp2*Mp
Cpq_p_3 = cp3*Mp
Cpq_c_1 = cc1*Mc
Cpq_c_2 = cc2*Mc
Cpq_c_3 = cc3*Mc
Cpq_a = ca*Ma

print('Cpq_p_1 =', Cpq_p_1)
print('Cpq_p_2 =', Cpq_p_2)
print('Cpq_p_3 =', Cpq_p_3)
print('Cpq_c_1 =', Cpq_c_1)
print('Cpq_c_2 =', Cpq_c_2)
print('Cpq_c_3 =', Cpq_c_3)
print('Cpq_a =', Cpq_a)

Cvq_p_1 = Cpq_p_1 - (9*(alphap**2)*kp*V0p*Tpm1)
Cvq_p_2 = Cpq_p_2 - (9*(alphap**2)*kp*V0p*Tpm2)
Cvq_p_3 = Cpq_p_3 - (9*(alphap**2)*kp*V0p*Tpm3)
Cvq_c_1 = Cpq_c_1 - (9*(alphac**2)*kc*V0c*Tcm1)
Cvq_c_2 = Cpq_c_2 - (9*(alphac**2)*kc*V0c*Tcm2)
Cvq_c_3 = Cpq_c_3 - (9*(alphac**2)*kc*V0c*Tcm3)
Cvq_a = Cpq_a - (9*(alphaa**2)*ka*V0a*Tam)

print('Cvq_p_1 =', Cvq_p_1)
print('Cvq_p_2 =', Cvq_p_2)
print('Cvq_p_3 =', Cvq_p_3)
print('Cvq_c_1 =', Cvq_c_1)
print('Cvq_c_2 =', Cvq_c_2)
print('Cvq_c_3 =', Cvq_c_3)
print('Cvq_a =', Cvq_a)

cp = 0.129
cc = 0.715
cal = 0.896

Cpl_p = cp*Mp
Cpl_c = cc*Mc
Cpl_a = cal*Ma

Cvl_p1 = Cpl_p - (9*(alphap**2)*kp*V0p*Tpm1)
Cvl_p2 = Cpl_p - (9*(alphap**2)*kp*V0p*Tpm2)
Cvl_p3 = Cpl_p - (9*(alphap**2)*kp*V0p*Tpm3)
Cvl_c1 = Cpl_c - (9*(alphac**2)*kc*V0c*Tcm1)
Cvl_c2 = Cpl_c - (9*(alphac**2)*kc*V0c*Tcm2)
Cvl_c3 = Cpl_c - (9*(alphac**2)*kc*V0c*Tcm3)
Cvl_a = Cpl_a - (9*(alphaa**2)*ka*V0a*Tam)

print('Mit Literaturwerten:')

print('Cpl_p =', Cpl_p)
print('Cpl_c =', Cpl_c)
print('Cpl_a =', Cpl_a)

print('Cvl_p1 =', Cvl_p1)
print('Cvl_p2 =', Cvl_p2)
print('Cvl_p3 =', Cvl_p3)
print('Cvl_c1 =', Cvl_c1)
print('Cvl_c2 =', Cvl_c2)
print('Cvl_c3 =', Cvl_c3)
print('Cvl_a =', Cvl_a)

print('3R =', 3*R)

#np.savetxt("ErgebnisseB.txt", np.column_stack([a.n, 1/a.n, Tx.n, Ty.n, Txy.n, Tpw1.n, Tpk1.n,
#Tpm1.n, Tpw2.n, Tpk2.n, Tpm2.n, Tpw3.n, Tpk3.n, Tpm3.n, Tcw1.n, Tck1.n, Tcm1.n, Tcw2.n, Tck2.n, Tcm2.n,
#Tcw3.n, Tck3.n, Tcm3.n, Taw.n, Tak.n, Tam.n, cgmg.n, cp1.n, cp2.n, cp3.n, cc1.n, cc2.n, cc3.n, ca.n, a.s, 1/a.s,
#Tx.s, Ty.s, Txy.s, Tpw1.s, Tpk1.s, Tpm1.s, Tpw2.s, Tpk2.s, Tpm2.s, Tpw3.s, Tpk3.s, Tpm3.s, Tcw1.s, Tck1.s,
#Tcm1.s, Tcw2.s, Tck2.s, Tcm2.s, Tcw3.s, Tck3.s, Tcm3.s, Taw.s, Tak.s, Tam.s, cgmg.s, cp1.s, cp2.s, cp3.s,
#cc1.s, cc2.s, cc3.s, ca.s]))
