# Stage-2 envelope trade-off numbers (code units: Rsun, km/s).
import numpy as np
G=6.674e-8; Ms=1.989e33; Rs=6.957e10
gu=G*Ms/(Rs*1e10)            # G Msun in (km/s)^2 Rsun
Ma,Md,a=6.24,1.69,30.3; R=4.06; rout=16.241
om=np.sqrt(gu*(Ma+Md)/a**3); gm=gu*Ma
cph=15.5; g=gm/R**2
print("Omega code",om,"P_orb",2*np.pi/om,"g",g,"GM/R",gm/R,"vK",np.sqrt(gm/R))
Hph=cph**2/g; print("H_p(photosphere) Rsun",Hph,"H/R",Hph/R)
dlog=np.log(rout/R)/448; drR=R*dlog; print("stage-1 dlnr",dlog,"dr(R)",drR,"dr/H",drR/Hph)
# stage-1 measured: dt 6.58e-6 code, 1.29e5 cyc/orbit, 23 min/orbit/node, 448x4x2048
for nenv in (116,):
  rin=R*np.exp(-nenv*dlog); print("n_env",nenv,"r_in",rin,"r_in/R",rin/R,"nx1",448+nenv)
for s in (1.0,7.2):
  for npol in (1.5,3.0):
    for f in (0.8,0.7,0.6):
      r=f*R; psi=gm*(1/r-1/R) - (s*s-1)*om**2*(r*r-R*R)/2
      c2=cph**2+psi/(npol+1); cad=np.sqrt(5/3*c2); vphi=(s-1)*om*r
      dr=r*dlog; dt=0.3*dr/(cad+vphi)
      rho=(c2/cph**2)**npol
      print(f"s {s} n {npol} r/R {f}: c_iso {np.sqrt(c2):.0f} c_ad {cad:.0f} vphi {vphi:.0f} dt {dt:.2e} rho/rho_ph {rho:.2e}")
# depth below which H_p >= 4 dr (n=3): H = (cph^2 + g d/(n+1))/g
for npol in (1.5,3.0):
  d=(4*drR - Hph)*(npol+1); print("n",npol,"depth where H_p=4dr:",d,"Rsun",d/R,"R")
# option A: resolved isothermal skin
drA=Hph/4; print("A: dr",drA,"dt ratio",drR/drA, "e-folds 0.7R->R", gm*(1/(0.7*R)-1/R)/cph**2)
