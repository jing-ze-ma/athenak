"""Convert L. Ma's profile (r T rho M kappa F_r, cgs) to the 7-column he_star_m1 IC layout that
bsg_figs.fig3 reads (r rho eint F_r E T fmlt), and extend it above rmax as he_ic_linhao does
(8e018272): T = T(rmax), rho = max(rho_rmax exp(-g (r - rmax) kT/(mu m_u)), dfloor), F_r = 0,
(bsgli_1104: gext invsq, rho = rho_rmax exp(-(GM/cs2)(1/rmax - 1/r)) instead of his (r/rmax)^2 g).  Unperturbed state (the run adds the 5 % rho perturbation).
usage: make_ic_linhao.py PROFILE OUT"""
import sys
import numpy as np
A_RAD, KB, MU, MU_U, G, DFLOOR = 7.5657e-15, 1.380649e-16, 0.62, 1.66053907e-24, 6.67408e-8, 1e-20
p = np.loadtxt(sys.argv[1])
r, T, rho, M, F = p[:, 0], p[:, 1], p[:, 2], p[:, 3], p[:, 5]
rmax, Tm, rm, Mm = r[-1], T[-1], rho[-1], M[-1]
ra = np.linspace(rmax, 5.6e12, 4001)[1:]          # to just above 80 Rsun
GM = 2.6542488e27   # linhao-prod-1104 he_ic_linhao_gext = invsq: hydrostatic in the run's GM/r^2
cs2 = KB * Tm / (MU * MU_U)
rha = np.maximum(rm * np.exp(-(GM / cs2) * (1.0 / rmax - 1.0 / ra)), DFLOOR)
R = np.concatenate([r, ra]); RHO = np.concatenate([rho, rha])
TT = np.concatenate([T, np.full_like(ra, Tm)]); FF = np.concatenate([F, np.zeros_like(ra)])
EINT = 1.5 * RHO * KB * TT / (MU * MU_U)
out = np.column_stack([R, RHO, EINT, FF, A_RAD * TT ** 4, TT, np.zeros_like(R)])
np.savetxt(sys.argv[2], out, fmt='%.10e', header=(
    'BSG page overlay IC: L. Ma linhao_profile_full.txt (he_ic_linhao, unperturbed) in he_star_m1 7-col '
    'layout; above rmax %.6e the he_ic_linhao isothermal extension\nr[cm] rho eint F_r E=aT^4 T[K] fmlt(0)' % rmax))
print('rows', len(R), 'rmax', rmax, 'T', Tm, 'rho', rm)
