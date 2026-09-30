"""Exact isothermal emission / heating reference for a ck_dump_kap column dump.
usage: python3 exact.py dump.txt [tag]
Columns (see two_stream_rt.hpp dump): 0 i 1 r_face 2 p 3 T 4 F_lw 10 Src 13 rho
14..24 F_b 25..35 B_b  41..128 kappa(b,g)  129..136 g weights.  Row i = face i / cell i."""
import sys
import numpy as np
NB, NG = 11, 8
d = np.loadtxt(sys.argv[1])
tag = sys.argv[2] if len(sys.argv) > 2 else ''
rf = d[:, 1]                       # faces is..ie+1
nc = len(rf) - 1
rho = d[:nc, 13]
T = d[:nc, 3]
Fb = d[:, 14:25]                   # per face
Bb = d[:nc, 25:36]
kap = d[:nc, 41:41+NB*NG].reshape(nc, NB, NG)
gw = d[0, 41+NB*NG:41+NB*NG+NG]
src = d[:nc, 10]
rin, rtop = rf[0], rf[-1]
kr = kap*rho[:, None, None]        # 1/cm


def seg(p, r1, r2):
    """path length inside each shell between radii r1<r2 on a ray of impact p."""
    a = np.clip(rf[:-1], r1, r2)
    b = np.clip(rf[1:], r1, r2)
    a = np.maximum(a, p)
    b = np.maximum(b, p)
    return np.sqrt(b*b - p*p) - np.sqrt(a*a - p*p)


# luminosity: impact-parameter quadrature, dense near every face
pts = np.unique(np.concatenate([np.linspace(0, rin, 200)] +
                [np.linspace(rf[i], rf[i+1], 60) for i in range(nc)]))
pm = 0.5*(pts[1:] + pts[:-1])
dp2 = pts[1:]**2 - pts[:-1]**2          # int 2 p dp over the bin
Lb = np.zeros(NB)
tau_top = np.zeros((NB, NG))
for p, w2 in zip(pm, dp2):
    if p < rin:
        Lb += 4*np.pi**2*Bb[0]*w2      # isothermal: hits the (isothermal) wall -> I = B
        continue
    s = 2*seg(p, p, rtop)
    tau = np.einsum('i,ibg->bg', s, kr)
    Lb += 4*np.pi**2*Bb[0]*w2*np.einsum('g,bg->b', gw, 1 - np.exp(-tau))
Lcode = 4*np.pi*rtop**2*Fb[-1]
Lph = 4*np.pi*rin**2*np.pi*Bb[0]
print(f'# {tag} T={T.mean():.1f} (spread {T.max()-T.min():.2e}) rin={rin:.5e} rtop={rtop:.5e}'
      f' x=A_top/A_in={(rtop/rin)**2:.4f}')
print('# sum pi B / sigma T^4 =', np.pi*Bb[0].sum()/(5.670374e-5*T[0]**4))
for b in range(NB):
    print(f'band {b:2d}  Lcode/Lexact {Lcode[b]/Lb[b]:.4f}   Lexact/(A_in pi B) {Lb[b]/Lph[b]:.4f}')
print(f'TOTAL Lcode/Lexact {Lcode.sum()/Lb.sum():.4f}  Lexact/(A_in pi B) {Lb.sum()/Lph.sum():.4f}'
      f'  Lcode/(A_in pi B) {Lcode.sum()/Lph.sum():.4f}')
# vertical optical depth from each face to the top, band-g min/median
tv = np.array([np.einsum('i,ibg->bg', np.diff(rf)*(np.arange(nc) >= i), kr) for i in range(nc)])
# exact mean intensity at cell centres -> net heating 4 pi kappa rho (J - B)
if len(sys.argv) > 3:
    rc = 0.5*(rf[1:] + rf[:-1])
    mu = np.polynomial.legendre.leggauss(48)
    icells = range(int(sys.argv[3]), nc)
    print('# i  r/rin  tau_v(min g) tau_v(max g)  Q_code  Q_exact  ratio   [erg/s/cm3]')
    for i in icells:
        r = rc[i]
        J = np.zeros((NB, NG))
        for m0, w0 in zip(*mu):
            for sgn in (1, -1):
                mm = 0.5*(m0 + 1)*sgn              # mu in (0,1) or (-1,0)
                ww = 0.5*w0*0.5                    # (1/2) int_-1^1 -> weights
                p = r*np.sqrt(1 - mm*mm)
                if sgn > 0:                        # outward: came from inside
                    if p < rin:
                        J += ww*1.0
                        continue
                    s = seg(p, p, r) + seg(p, p, rtop)
                else:                              # inward: came from the top
                    s = seg(p, r, rtop)
                tau = np.einsum('i,ibg->bg', s, kr)
                J += ww*(1 - np.exp(-tau))
        Qex = np.sum(4*np.pi*kr[i]*Bb[i][:, None]*(J - 1)*gw[None, :])
        print(f'{i:3d} {r/rin:.4f} {tv[i].min():.3e} {tv[i].max():.3e} {src[i]: .4e} {Qex: .4e}'
              f' {src[i]/Qex if Qex != 0 else np.nan: .4f}')
if rho.min() < 1e-10:
    iph = int(np.argmax(rho < 1e-10))
    rph = rf[iph]
    xx = (rtop/rph)**2
    Aph = 4*np.pi*rph**2*np.pi*Bb[0].sum()
    print(f'SHARP r_ph={rph:.5e} x={xx:.4f} Lexact/(A_ph pi B)={Lb.sum()/Aph:.5f}'
          f' Lcode/(A_ph pi B)={Lcode.sum()/Aph:.5f} pred 2x/(1+x)={2*xx/(1+xx):.5f}'
          f' tau_v shell max={tv[iph].max():.2e}  tau_v(ph cell, min g)={(kr[iph-1]*(rf[iph]-rf[iph-1])).min() if iph>0 else 0:.2e}')
