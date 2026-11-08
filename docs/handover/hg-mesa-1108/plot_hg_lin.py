#!/usr/bin/env python3
"""hg-mesa-1108: plot the HG inner-boundary luminosity history from hg_lin.csv.

CSV columns (header row, one row per MESA model from the profile14 match onward):
  model, star_age_yr, t_cc_yr, t_si_yr, L_surf, R_surf, Teff,
  L_r3, L_r8, L_m3, r_m3, L_m8, r_m8, L_He3a, r_He3a, L_Hebase, r_Hebase, log_Tc
L in Lsun, r and R in Rsun, Teff and Tc in K; t_cc_yr = age(final) - age,
t_si_yr = age(Si ignition) - age (negative after ignition).
Usage: python3 plot_hg_lin.py hg_lin.csv [out.png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

f = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "hg_lin.png"
d = np.genfromtxt(f, delimiter=",", names=True)
t = np.maximum(d["t_cc_yr"], 1e-6)

fig, ax = plt.subplots(3, 1, figsize=(7, 10), sharex=True)
a = ax[0]
for key, lab in [("L_r8", "L(r = 8 Rsun) [deck r_in]"), ("L_r3", "L(r = 3 Rsun)"),
                 ("L_m8", "L(m = m(8 Rsun, p14))"), ("L_m3", "L(m = m(3 Rsun, p14))"),
                 ("L_surf", "L surface"), ("L_He3a", "L at max eps_3a"),
                 ("L_Hebase", "L at He-envelope base")]:
    a.plot(t, d[key], label=lab)
a.set_yscale("log")
a.set_ylabel("L [Lsun]")
a.legend(fontsize=7)
a = ax[1]
for key, lab in [("R_surf", "R surface"), ("r_m8", "r(m_in8)"), ("r_m3", "r(m_in3)"),
                 ("r_He3a", "r of max eps_3a"), ("r_Hebase", "r He-envelope base")]:
    a.plot(t, d[key], label=lab)
a.set_yscale("log")
a.set_ylabel("r [Rsun]")
a.legend(fontsize=7)
a = ax[2]
a.plot(t, d["Teff"], label="Teff")
a.set_ylabel("Teff [K]")
a.set_xscale("log")
a.invert_xaxis()
a.set_xlabel("time to core collapse [yr]")
ts = d["t_si_yr"]
if np.any(ts <= 0):
    tsi = t[np.argmax(ts <= 0)]
    for a in ax:
        a.axvline(tsi, color="k", ls=":", lw=0.8)
fig.tight_layout()
fig.savefig(out, dpi=130)
print("wrote", out)
