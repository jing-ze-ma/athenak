"""system.py (dhj_planets_1101/wasp189b; copied from toi1518b): adopted literature parameters -> derived numbers
(Teq, omega, R_p, a, Roche sizes) and the tide/gravity ratio at given radii.
usage: system.py [r1 r2 ...]   (radii [cm] for the tide table, e.g. ap R_p x1max)"""
import sys
import numpy as np

G = 6.67430e-8
GMJ = 1.26686534e23          # IAU 2015 nominal (planet_setup.py MJ)
GMSUN = 1.3271244e26
RJ, RSUN, AU, DAY = 7.1492e9, 6.957e10, 1.495978707e13, 86400.0
# Lendl et al. 2020, A&A 643, A94 (CHEOPS + HARPS/CORALIE RV; NASA archive default)
P = 2.7240330                # d
ARS = 4.600                  # a/R_*
RPRS = 0.07045               # R_p/R_* (CHEOPS)
RS = 2.362                   # R_sun
MS = 2.030                   # M_sun
TEFF = 8000.0
MP = 1.99                    # M_J
A_AU = 0.05053
RP_RJ = 1.619                # Lendl Table: R_p (adopted for --rp)

rs = RS*RSUN
a = ARS*rs
rp = RPRS*rs
gmp = MP*GMJ
gms = MS*GMSUN
q = gmp/gms
om = 2*np.pi/(P*DAY)
teq = TEFF*np.sqrt(1.0/(2.0*ARS))
print('R_p = RPRS R_* = %.6e cm = %.4f R_J (Lendl+2020: 1.619 +- 0.021)' % (rp, rp/RJ))
print('a = ARS R_* = %.6e cm = %.5f AU (Lendl a: %.5f AU)' % (a, a/AU, A_AU))
print('Kepler a from P, M_*+M_p: %.5f AU' % ((G*0 + (gms + gmp)*(P*DAY/2/np.pi)**2)**(1/3)/AU))
print('Teq = Teff sqrt(R_*/2a) = %.1f K (zero albedo, full redistribution)' % teq)
print('omega = 2 pi / P = %.6e s^-1 (P = %.1f s)' % (om, P*DAY))
print('g(R_p) = %.1f cm/s^2, log g = %.3f' % (gmp/rp**2, np.log10(gmp/rp**2)))
print('q = M_p/M_* = %.4e' % q)
rh = a*(q/3)**(1/3)
q3, q23 = q**(1/3), q**(2/3)
rl = a*0.49*q23/(0.6*q23 + np.log(1 + q3))
print('Hill (~L1) radius a (q/3)^(1/3) = %.4e cm = %.3f R_p' % (rh, rh/rp))
print('Eggleton Roche-lobe volume radius = %.4e cm = %.3f R_p' % (rl, rl/rp))
print('Stellar rotation: v sin i 93.1 km/s (Anderson+2018) -> P_rot,* <= %.2f d'
      % (2*np.pi*rs/93.1e5/DAY))
radii = [float(v) for v in sys.argv[1:]]
if radii:
    print('\n tide/gravity ratios (point mass g = GM_p/r^2):')
    print('   r [cm]      r/R_p   r/R_L1   2(r/a)^3/q (stellar tide, substellar)  '
          '3(r/a)^3/q (tide + orbital centrifugal, substellar)  om^2 r/g (spin '
          'centrifugal, equator; in rot_potential)')
    for r in radii:
        g = gmp/r**2
        print('  %.4e  %6.3f  %6.3f   %.3e   %.3e   %.3e'
              % (r, r/rp, r/rh, 2*(r/a)**3/q, 3*(r/a)**3/q, om**2*r/g))
