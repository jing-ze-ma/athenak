"""system.py (dhj_planets_1101): adopted literature parameters -> derived numbers (Teq, omega,
g, a, Roche/Hill sizes), tide/gravity ratios at given radii, dynamo estimate.
usage: system.py PLANET [r1 r2 ...]   (PLANET = wasp76b | wasp178b; radii [cm], e.g. ap x1max)
Same formulas as toi1518b/scripts/system.py; parameters in the SYS table (README section 1)."""
import sys
import numpy as np

G = 6.67430e-8
GMJ = 1.26686534e23          # IAU 2015 nominal (planet_setup.py MJ)
GMSUN = 1.3271244e26
RJ, RSUN, AU, DAY = 7.1492e9, 6.957e10, 1.495978707e13, 86400.0
SIG = 5.670374419e-5
SYS = {
    # Ehrenreich+2020 Nature 580, 597 (ED Table; NASA archive default)
    'wasp76b': dict(P=1.80988198, ARS=4.08, RP=1.854, RS=1.756, MS=1.458, TEFF=6329.0,
                    MP=0.894, A_AU=0.0330, TEQ_PAPER=2229.0, VSINI=1.48),
    # Lothringer+2025 AJ 169, 274 (new stellar characterisation + joint HST/JWST/TESS fit:
    # Teff, R*, M*, a/R*, Rp); P Damasceno+2024 A&A 689, A54; Mp Hellier+2019 (RV)
    'wasp178b': dict(P=3.3448371, ARS=7.037, RP=1.77, RS=1.70, MS=2.02, TEFF=9200.0,
                     MP=1.66, A_AU=0.0561, TEQ_PAPER=2398.0, VSINI=10.48),
}
name = sys.argv[1]
s = SYS[name]
rs = s['RS']*RSUN
a = s['ARS']*rs
rp = s['RP']*RJ
gmp = s['MP']*GMJ
gms = s['MS']*GMSUN
q = gmp/gms
P = s['P']*DAY
om = 2*np.pi/P
teq = s['TEFF']*np.sqrt(1.0/(2.0*s['ARS']))
print('%s' % name)
print('R_p = %.6e cm = %.4f R_J;  R_p/R_* = %.5f' % (rp, s['RP'], rp/rs))
print('a = ARS R_* = %.6e cm = %.6f AU (quoted a: %.4f AU); Kepler a(P, M_*+M_p) = %.5f AU'
      % (a, a/AU, s['A_AU'], ((gms + gmp)*(P/2/np.pi)**2)**(1/3)/AU))
print('Teq = Teff sqrt(R_*/2a) = %.1f K (zero albedo, full redistribution); Seidel+2026 '
      'ED Table 1: %.0f K (%+.1f %%)' % (teq, s['TEQ_PAPER'], 100*(teq/s['TEQ_PAPER'] - 1)))
F = 4*SIG*teq**4
tint = max(0.39*teq*np.exp(-(np.log10(F) - 9.0 - 0.14)**2/1.095), 100.0)
print('T_int (Thorngren+2019 fit as rce_w121.py and two_stream_rt::get_Tint) = %.1f K' % tint)
print('omega = 2 pi / P = %.6e s^-1 (P = %.1f s); v_eq(R_p) = %.3f km/s'
      % (om, P, om*rp/1e5))
print('g(R_p) = %.1f cm/s^2, log g = %.3f' % (gmp/rp**2, np.log10(gmp/rp**2)))
print('q = M_p/M_* = %.4e' % q)
rh = a*(q/3)**(1/3)
q3, q23 = q**(1/3), q**(2/3)
rl = a*0.49*q23/(0.6*q23 + np.log(1 + q3))
print('Hill (~L1) radius = %.4e cm = %.3f R_p;  Eggleton Roche-lobe radius = %.4e cm = %.3f R_p'
      % (rh, rh/rp, rl, rl/rp))
print('Stellar rotation: v sin i %.2f km/s -> P_rot,* sin^-1 i <= %.2f d'
      % (s['VSINI'], 2*np.pi*rs/(s['VSINI']*1e5)/DAY))
# dynamo: B ~ F_int^(1/3) relative to Jupiter (Christensen+2009 energy-flux scaling, density
# factor and dynamo-region geometry ignored).  Jupiter: F_int = 7.485 W/m^2 (Li+2018, T_int
# 107 K); older 5.44 W/m^2 (Hanel+1981, 99 K). B_J: dipole equatorial surface 4.2 G
# (rms surface ~ 7 G, Juno max ~ 20 G, Connerney+2018/2022).
for fj, lab in ((7.485e3, 'Li+2018'), (5.44e3, 'Hanel+1981')):
    r = (SIG*tint**4/fj)**(1/3)
    print('dynamo: F_int/F_J(%s) = %.0f -> B/B_J = %.2f -> B_dip,eq ~ %.0f G (B_J 4.2 G), '
          'rms ~ %.0f G (7 G)' % (lab, SIG*tint**4/fj, r, 4.2*r, 7*r))
radii = [float(v) for v in sys.argv[2:]]
if radii:
    print('\n tide/gravity ratios (point mass g = GM_p/r^2):')
    print('   r [cm]      r/R_p   r/R_L1   2(r/a)^3/q (stellar tide, substellar)  '
          '3(r/a)^3/q (tide + orbital centrifugal)  om^2 r/g (spin centrifugal, equator; '
          'in rot_potential)')
    for r in radii:
        g = gmp/r**2
        print('  %.4e  %6.3f  %6.3f   %.3e   %.3e   %.3e'
              % (r, r/rp, r/rh, 2*(r/a)**3/q, 3*(r/a)**3/q, om**2*r/g))
