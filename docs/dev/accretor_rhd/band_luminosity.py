# S3a: the luminosity of the equatorial band (2-D slab) of the Roche-distorted gainer, as a
# fraction of the column's L = 4 pi r90^2 F90, when F ~ g on each equipotential (diffusion,
# von Zeipel) and T = T(psi): ratio = <(g/g90) (r/r90)^2 / cos(alpha)>_phi.  Explains
# Lmeas/L = 0.964 of the static radiative run (prediction 0.9605 at the photosphere).
import numpy as np
from scipy.optimize import brentq
kG,kMsun,kRsun,kVel=6.674e-8,1.989e33,6.957e10,1e5
gu=kG*kMsun/(kRsun*kVel**2); ma,md,a=16.,18.,32.5575
gma,gmd=gu*ma,gu*md; om=np.sqrt(gu*(ma+md)/a**3); xcm=a*md/(ma+md)
def P(r,ph):
    c=np.cos(ph); rd=np.sqrt(r*r+a*a-2*a*r*c)
    return -gma/r-gmd/rd-0.5*om*om*(r*r+xcm*xcm-2*xcm*r*c)
racc=9.00129260997162
for psi in (0.0, 2.0e4, 6.0e4, 1.4343e5):
    ph0=P(racc,np.pi/2)+psi* -1  # Phi = Phi_s - psi
    target=P(racc,np.pi/2)-psi
    phs=np.linspace(0,2*np.pi,4001)[:-1]
    r=np.array([brentq(lambda x:P(x,p)-target,3,13.4) for p in phs])
    h=1e-6
    dr=np.array([(P(x+h,p)-P(x-h,p))/(2*h) for x,p in zip(r,phs)])
    dp=np.array([(P(x,p+h)-P(x,p-h))/(2*h)/x for x,p in zip(r,phs)])
    g=np.sqrt(dr**2+dp**2); cosa=dr/g
    k=np.argmin(abs(phs-np.pi/2)); g90=g[k]; r90=r[k]
    # flux through the band surface per dtheta dphi: F_n r^2/cos(a), F_n = F90 g/g90
    ratio=np.mean((g/g90)*(r/r90)**2/cosa)
    print('psi %.4g r90 %.4f r(phi0) %.4f  L_band/L(column) = %.5f' % (psi, r90, r[0], ratio))
