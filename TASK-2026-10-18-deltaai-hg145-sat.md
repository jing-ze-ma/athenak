# TASK raven -> DeltaAI (2026-10-18): He giant -- saturation check of the deep convection on your 0-28 d dumps (CPU only, ~15 min)

Raven diagnosis: the deep KE growth is the Fe-bump convection zone (base ~19 Rsun, unstable in our model and in MESA, N^2 ~ -5e-12)
turning on from the 1e-3 seed; e-fold of KE 1.5 -> 2.9 d and slowing; MLT saturation estimate KE ~1-2e42 erg (wedge), v ~10-12 km/s.
Please confirm or refute with your dumps (reduce.py already writes Fenth/Fke per dump in data/red_*.npz):
1. Per dump at 2-d steps (and the last 3 dumps): L_conv,res(r) = Omega r^2 (Fenth + Fke) shell mean, L_rad(r) from the m1 flux,
   and the ratio L_conv,res / (L_in - L_rad) at 20-30, 30-60, 60-100 Rsun. Also KE_int e-folding time per 2-d window (hst col 11).
2. rms v_r and max |v_r| per band 8-19, 19-30, 30-60, 60-100 Rsun at the same times.
Verdict wanted: is L_conv,res/(L_in - L_rad) approaching 1 and the KE growth rate falling (= saturating)? Flag if KE > ~3e42 erg or
v_rms > ~20 km/s.
Report: NOTE-2026-10-18-deltaai-hg145-sat.md (table + 3-line verdict). No GPU, no new run.
